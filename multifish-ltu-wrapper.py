import sys,os,glob
import memoryCC as mcc


# ===================================================================================
# Global Vars
# ===================================================================================
# log is actually never set by the obj C side.
global log
log = True

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

def BlankLTUParameterDict():
    # Blank dict of LTU parameters that we use for initialisation.
    # Note that Daniel previously did dict(LTUParameterDict) where LTUParameterDict was a variable.
    # This only did a shallow copy, which meant the list entries in the blank dict were being reused across fish!
    return { 'resampledSequences' : [],
              'periodHistory' : [],
              'driftHistory' : [],
              'shifts' : [],
              'knownPhaseIndex' : -1,
              'knownPhase' : 0 }

# The nested dict / oracle containing the LTU parameters for multiple fish.
# At startup we will not have any entries, but at least one should be added by the Spim GUI during its own startup
multifishOracle = dict()

# helper functions that are not called by the LTU App
def isFishProfileInOracle(uniqueFishID):
    return (uniqueFishID in multifishOracle.keys())

def addFishToOracle(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == False):
        multifishOracle[uniqueFishID] = BlankLTUParameterDict()
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} is already in oracle')
    sys.stdout.flush()
    
def addFishToOracleIfNeeded(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == False):
        print(f'Adding unique fish ID {uniqueFishID} to the oracle')
        multifishOracle[uniqueFishID] = BlankLTUParameterDict()
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

def referencePhaseWasActivelySetForMostRecentSequence(fractionThroughSequence, uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        print(f'Updating knownPhase for unique fish ID {uniqueFishID}')
        multifishOracle[uniqueFishID]['knownPhaseIndex'] = len(multifishOracle[uniqueFishID]['resampledSequences']) - 1
        multifishOracle[uniqueFishID]['knownPhase'] = fractionThroughSequence * numSamplesPerPeriod
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} is not in oracle. Will add new entry with null parameters')
        addFishToOracle(uniqueFishID)
    sys.stdout.flush()

def updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        print(f'Updating LTU parameters for unique fish ID {uniqueFishID}')
        multifishOracle[uniqueFishID]['resampledSequences'] = resampledSequences
        multifishOracle[uniqueFishID]['periodHistory'] = periodHistory
        multifishOracle[uniqueFishID]['driftHistory'] = driftHistory
        multifishOracle[uniqueFishID]['shifts'] = shifts
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} is not in oracle. Will add new entry and update parameters')
        addFishToOracle(uniqueFishID)
        updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID)
    sys.stdout.flush()

def get4LTUParameters(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        ltuTuple = tuple(multifishOracle[uniqueFishID][key] for key in ['resampledSequences', 'periodHistory','driftHistory', 'shifts'])
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} was queried but is not in oracle. Will add new entry and return requested (empty) parameters')
        addFishToOracle(uniqueFishID)
        ltuTuple = get4LTUParameters(uniqueFishID)
    sys.stdout.flush()
    return ltuTuple

def get6LTUParameters(uniqueFishID):
    if (isFishProfileInOracle(uniqueFishID) == True):
        ltuTuple = tuple(multifishOracle[uniqueFishID][key] for key in ['resampledSequences', 'periodHistory','driftHistory', 'shifts', 'knownPhaseIndex', 'knownPhase'])
    else:
        print(f'Unexpected: unique fish ID {uniqueFishID} was queried but is not in oracle. Will add new entry and return requested (empty) parameters')
        addFishToOracle(uniqueFishID)
        ltuTuple = get6LTUParameters(uniqueFishID)
    sys.stdout.flush()
    return ltuTuple


# ===================================================================================
# Wrapper functions for the MemoryCC module
# ===================================================================================
'''
These functions call their respective functions in memoryCC for aligning reference sequences for timelapse imaging with 
the addition of interfacing with the multifish oracle. Each of these functions follows a common pattern on being called: 
- get the LTUparameters from the oracle
- combine with new data coming in from the Obj C side and pass to the old MemoryCC functions
- update LTUparameters in the oracle
- return required parameters back to the Obj C side.
'''

def getFractionalPhaseByAligningReferenceSequence(rawFrames, thisPeriod, thisDrift, maxOffsetToConsider, uniqueFishID):
    print(f'getFractionalPhaseByAligningReferenceSequence for unique fish ID {uniqueFishID}')
    ltuParameters = get6LTUParameters(uniqueFishID)
    resampledSequences, periodHistory, driftHistory, shifts, shiftSolution, _ = mcc.processNewReferenceSequence(rawFrames, thisPeriod, thisDrift, *ltuParameters, numSamplesPerPeriod, maxOffsetToConsider)
    print(f'getFractionalPhaseByAligningReferenceSequence completed for unique fish ID {uniqueFishID} (result {shiftSolution:.3f}, frac {(shiftSolution/numSamplesPerPeriod)%1.0:.3f})')
    updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID)
    _ltuParameters = get6LTUParameters(uniqueFishID)
    # Note that we never actually use the residuals that get returned.
    # Only the shiftSolution is actually need by the LTU helper app
    sys.stdout.flush()
    return (shiftSolution / numSamplesPerPeriod) % 1.0

def trimLTUHistory(trimToLength, uniqueFishID):
    print(f'Trim LTU history for fish {uniqueFishID}')
    ltuParameters = get4LTUParameters(uniqueFishID)
    returnTuple = mcc.trimLTUHistory(*ltuParameters, trimToLength)
    updateLTUParameters(*returnTuple, uniqueFishID)
    sys.stdout.flush()

def RoIForReferenceHistory(uniqueFishID):
    # this function only needs a reference to resampledSequences so I'm not going to call the updater
    # The tuple (-1,-1) is returned if the length of resampledSequences we pass in is zero.
    # If the fish ID doesn't exist then I think we should still create an entry in the oracle then recall the function.
    # This will still return (-1,-1) back to the obj C side BUT we wont crash by reading a non existent entry in the oracle.
    ltuParameters = get4LTUParameters(uniqueFishID)
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
    resampledSequences,_,_,_ = get4LTUParameters(uniqueFishID)
    return len(resampledSequences)

def resetRefFrameHistory(uniqueFishID):
    # Clear the parameters for this fish's entry in the oracle
    print(f'Reset ref frame history for fish {uniqueFishID} (was {numRefFrameSetsInHistory(uniqueFishID)} in history)')
    multifishOracle[uniqueFishID] = BlankLTUParameterDict()
    sys.stdout.flush()
