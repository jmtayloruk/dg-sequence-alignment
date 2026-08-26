import sys,os,glob
import memoryCC as mcc


# ===================================================================================
# Global Vars
# ===================================================================================
# log is actually never set by the obj C side.
global log
log = True

# numSamplesPerPeriod is hard coded in SyncControllerPythonInterface and never updated
global numSamplesPerPeriod
numSamplesPerPeriod = 60

# ===================================================================================
# The Multifish Oracle and Helper Functions
# ===================================================================================
'''
The Multifish oracle is just a nested dict where the inner dict tracks parameters of the long-term updater for a single fish. 
The outer dict collects these into the oracle^TM. For each fish theres is a key "N" where N is the uniqueFishID in spimGUI 
and the value is the inner dict of LTU parameters.
The Spim GUI should be responsible for making sure addFishToOracleIfNeeded is called whenever a new fish profile is created.
That means we can print an "Unexpected: ..." message if we ever find a uniqueFishID is missing from (or present in) the oracle
when we do not expect that. But we will handle that, so the Spim GUI shouldn't see any impact
(beyond potentially degradation of LTU performance for that fish).
Because the Spim GUI always ensures it has at least one fish profile, this oracle should operate transparently even when capturing a single-fish timelapse.
'''

# blank dict of LTU parameters that we only ever copy from and never update.
blankLTUParameterDict = { 'resampledSequences' : [],
                          'periodHistory' : [],
                          'driftHistory' : [],
                          'shifts' : []}

# The nested dict / oracle containing the LTU parameters for multiple fish.
# At startup we will not have any entries, but at least one should be added by the Spim GUI during its own startup
multifishOracle = dict()

# helper functions that are not called by the LTU App
def isFishProfileInOracle(uniqueFishID):
    return (uniqueFishID in multifishOracle.keys())

def addFishToOracle(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == False):
        multifishOracle[uniqueFishID] = dict(blankLTUParameterDict)
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} is already in oracle')
    sys.stdout.flush()
    
def addFishToOracleIfNeeded(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == False):
        print(f'Adding unique fish ID {uniqueFishID} to the oracle')
        multifishOracle[uniqueFishID] = dict(blankLTUParameterDict)
    else:
        print(f'For information: unique fish ID {uniqueFishID} is already in oracle')
    sys.stdout.flush()

def removeFishFromOracle(uniqueFishID):
    # To remove a fish from the oracle we remove the entry for that key.
    # Because the uniqueFishID is decoupled from the fishIndex in the Spim GUI,
    # we don't need to shuffle anything else around.
    uniqueFishIDs = sorted(keys for keys in multifishOracle.keys())
    numFishInOracle = len(uniqueFishIDs)
    if (isFishProfileInOracle(uniqueFishID) == False):
        # This fish ID is not in the oracle. That is not necessarily a problem,
        # it may just mean that no sync has been performed for that fish, even though
        # the fish was defined within the Spim GUI.
        print(f'Asked to delete unique fish ID {uniqueFishID} but it is not in oracle')
    else:
        # Delete this entry from the oracle
        print(f'Deleting unique fish ID {uniqueFishID} from oracle')
        del multifishOracle[uniqueFishID]
    sys.stdout.flush()

def updateLTUParameters(resampledSequences, periodHistory, driftHistory,  shifts, uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        print(f'Updating LTU parameters for unique fish ID {uniqueFishID}')
        parameterDict = {   'resampledSequences' : resampledSequences,
                            'periodHistory' : periodHistory,
                            'driftHistory' : driftHistory,
                            'shifts' : shifts
                         }
        multifishOracle[uniqueFishID] = parameterDict
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} is not in oracle. Will add new entry and update parameters')
        addFishToOracle(uniqueFishID)
        updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID)
    sys.stdout.flush()

def getLTUParameters(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        ltuTuple = tuple(multifishOracle[uniqueFishID][key] for key in ['resampledSequences', 'periodHistory','driftHistory', 'shifts'])
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} was queried but is not in oracle. Will add new entry and return requested (empty) parameters')
        addFishToOracle(uniqueFishID)
        ltuTuple = getLTUParameters(uniqueFishID)
    sys.stdout.flush()
    return ltuTuple


# ===================================================================================
# Wrapper functions for the MemoryCC module
# ===================================================================================
'''
These functions call their respective functions in memoryCC for aligning reference sequences for timelapse imaging with 
the addition of interfacing with the multifish oracle. Each of these functions follows a common pattern on being called: get the LTUparameters from the oracle;
combine with new data coming in from the Obj C side and pass to the old MemoryCC functions; update LTUparameters in the oracle; return required parameters back to the Obj C side.
'''

def processNewReferenceSequence(rawFrames, thisPeriod, thisDrift, knownPhaseIndex, knownPhase, maxOffsetToConsider, uniqueFishID):
    print(f'processNewReferenceSequence for unique fish ID {uniqueFishID}')
    ltuParameters = getLTUParameters(uniqueFishID)
    # we never actually use the residuals that get returned. Only the shiftSolution actually need by the LTU helper app
    resampledSequences, periodHistory, driftHistory, shifts, shiftSolution, _ = mcc.processNewReferenceSequence(rawFrames, thisPeriod, thisDrift, *ltuParameters, knownPhaseIndex, knownPhase, numSamplesPerPeriod, maxOffsetToConsider)
    updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID)
    print(f'processNewReferenceSequence completed for unique fish ID {uniqueFishID} (result {shiftSolution})')
    sys.stdout.flush()
    return shiftSolution

def trimLTUHistory(trimToLength, uniqueFishID):
    print(f'Trim LTU history for fish {uniqueFishID}')
    ltuParameters = getLTUParameters(uniqueFishID)
    returnTuple = mcc.trimLTUHistory(*ltuParameters, trimToLength)
    updateLTUParameters(*returnTuple,uniqueFishID)
    sys.stdout.flush()

def RoIForReferenceHistory(uniqueFishID):
    # this function only needs a reference to resampledSequences so I'm not going to call the updater
    # The tuple (-1,-1) is returned if the length of resampledSequences we pass in is zero.
    # If the fish ID doesn't exist then I think we should still create an entry in the oracle then recall the function.
    # This will still return (-1,-1) back to the obj C side BUT we wont crash by reading a non existent entry in the oracle.
    ltuParameters = getLTUParameters(uniqueFishID)
    roi = mcc.RoIForReferenceHistory(ltuParameters[0])
    sys.stdout.flush()
    return roi


# ===================================================================================
# Additional functions used by the LTU app
# ===================================================================================
'''
numRefFrameSetsInHistory and resetRefFrameHistory, were previously methods belonging to PythonService.
But since resampledSequences, periodHistory, driftHistory, and shifts are now only tracked on the python side,
resampledSequences doesn't exist in the Obj C.
So these functions are implemented here and are called by the Obj C side. 
'''

def numRefFrameSetsInHistory(uniqueFishID):
    # Return number of sets of reference sequences in the list of resampledSequences
    resampledSequences ,_ ,_, _= getLTUParameters(uniqueFishID)
    return len(resampledSequences)

def resetRefFrameHistory(uniqueFishID):
    # set the parameters for that entry in the oracle to an empty list
    print(f'Reset ref frame history for fish {uniqueFishID}')
    updateLTUParameters([],[],[],[], uniqueFishID)
