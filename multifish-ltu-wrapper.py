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
The outer dict collects these into the oracle^TM. For each fish theres is a key "N" where N is the fishIndex in spimGUI 
and the value is the inner dict of LTU parameters.
The various helper functions should ensure that behaviour is "mostly" consistent with how the settingsForFish container is handled
in the main spimGUI in that entries can be added and removed (controlled by the main app) BUT there is always an entry for 
fishIndex = 0 ("0"). This means you can still record a timelapse for a single fish.
'''

# blank dict of LTU parameters that we only ever copy from and never update.
LTUParameterDict = { 'resampledSequences' : [],
                     'periodHistory' : [],
                     'driftHistory' : [],
                     'shifts' : []}

# the nested dict / oracle containing the LTU parameters for multiple fish. There will always be an entry with key "0"
multifishOracle = {0 : dict(LTUParameterDict)}

# helper functions that are not called by the LTU App
def isFishProfileInOracle(fishIndex):
    return (fishIndex in multifishOracle.keys())

def addFishToOracle(fishIndex):
    if (isFishProfileInOracle(fishIndex) == False):
        multifishOracle[fishIndex] = dict(LTUParameterDict)

def removeFishFromOracle(fishIndex):
    # removing a fish from the oracle removes the entry at that key but also decrements the key number by one to match the behaviour of the spimGUI obj C side
    # im mostly just going to copy the logic of the SpimApplication.removeFish method
    fishIndices = sorted(keys for keys in multifishOracle.keys())
    maxFishIndex = max(fishIndices)
    numFishInOracle = len(fishIndices)
    if (isFishProfileInOracle(fishIndex) == False):
        # This fish index is not in the oracle. That is not necessarily a problem,
        # it may just mean that no sync has been performed for that fish, even though
        # the fish was defined within the Spim GUI.
        print(f"Fish index {fishIndex} is not in oracle.")
    elif (numFishInOracle == 1):
        # there is only one fish in the oracle. If all rules have been followed this will be index 0
        print(f"We can't delete the only entry in the oracle. Resetting the LTU parameters instead for fish index {fishIndex}")
        assert(fishIndex == 0)
        updateLTUParameters([],[],[],[], fishIndex)
    elif (fishIndex == maxFishIndex):
        # if its the final entry we want to remove just delete it
        del multifishOracle[fishIndex]
    else:
        # shuffle all the entries down by updating each entry with the LTU parameters from the entry above and deleting the final entry
        # again we're relying on nothing fishy happening and that there are continuous indices from fishIndex to maxIndex
        for k in range(fishIndex, maxFishIndex,1):
            updateLTUParameters(*getLTUParameters(k+1), k)
        # delete final entry because there is no successor to update LTUparameters from.
        del multifishOracle[maxFishIndex]

def updateLTUParameters(resampledSequences, periodHistory, driftHistory,  shifts, fishIndex):
    if(isFishProfileInOracle(fishIndex) == True):
        parameterDict = {   'resampledSequences' : resampledSequences,
                            'periodHistory' : periodHistory,
                            'driftHistory' : driftHistory,
                            'shifts' : shifts
                         }
        multifishOracle[fishIndex] = parameterDict
    else:
        print(f'Fish Index {fishIndex} is not in oracle. Will add new entry and update parameters')
        addFishToOracle(fishIndex)
        updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, fishIndex)

def getLTUParameters(fishIndex):
    if (isFishProfileInOracle(fishIndex) == True):
        ltuTuple = tuple(multifishOracle[fishIndex][key] for key in ['resampledSequences', 'periodHistory','driftHistory', 'shifts'])
    else:
        print(f'Fish Index {fishIndex} is not in oracle. Will add new entry and return requested parameters')
        addFishToOracle(fishIndex)
        ltuTuple = getLTUParameters(fishIndex)
    return ltuTuple


# ===================================================================================
# Wrapper functions for the MemoryCC module
# ===================================================================================
'''
These functions call their respective functions in memoryCC for aligning reference sequences for timelapse imaging with 
the addition of interfacing with the multifish oracle. Each of these functions follows a common pattern on being called: get the LTUparameters from the oracle;
combine with new data coming in from the Obj C side and pass to the old MemoryCC functions; update LTUparameters in the oracle; return required parameters back to the Obj C side.
'''

def processNewReferenceSequence(rawFrames, thisPeriod, thisDrift,  knownPhaseIndex,knownPhase, maxOffsetToConsider, fishIndex = 0):
    
    ltuParameters = getLTUParameters(fishIndex)
    # we never actually use the residuals that get returned. Only the shiftSolution actually need by the LTU helper app
    resampledSequences, periodHistory, driftHistory, shifts, shiftSolution, _ = mcc.processNewReferenceSequence(rawFrames, thisPeriod, thisDrift, *ltuParameters, knownPhaseIndex, knownPhase, numSamplesPerPeriod, maxOffsetToConsider)
    updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, fishIndex)
    return shiftSolution

def trimLTUHistory(trimToLength, fishIndex = 0):
    ltuParameters = getLTUParameters(fishIndex)
    returnTuple = mcc.trimLTUHistory(*ltuParameters, trimToLength)
    updateLTUParameters(*returnTuple,fishIndex)

def RoIForReferenceHistory(fishIndex = 0):
    # this function only needs a reference to resampledSequences so I'm not going to call the updater
    # The tuple (-1,-1) is returned if the length of resampledSequences we pass in is zero.
    # If the fish index doesn't exist then I think we should still create an entry in the oracle then recall the function.
    # This will still return (-1,-1) back to the obj C side BUT we wont crash by reading a non existent entry in the oracle.
    ltuParameters = getLTUParameters(fishIndex)
    roi = mcc.RoIForReferenceHistory(ltuParameters[0])
    return roi


# ===================================================================================
# Additional functions used by the LTU app
# ===================================================================================
'''
These are additional functions that either needed to be implemented for multifish imaging (resetMultifishOracle, and clearMultifishOracle)
The other two functions, numRefFrameSetsInHistory and resetRefFrameHistory, were previously methods belonging to PythonService.
But since resampledSequences, periodHistory, driftHistory, and shifts are now only tracked on the python side then resampledSequences doesn't exist in the Obj C.
So these functions are implemented here and are called by the Obj C side. The Obj C side requires each function return something so I return 1. 
'''

def numRefFrameSetsInHistory(fishIndex = 0):
    # Return number of sets of reference sequences in the list of resampledSequences
    resampledSequences ,_ ,_, _= getLTUParameters(fishIndex)
    return len(resampledSequences)

def resetRefFrameHistory(fishIndex = 0):
    # set the parameters for that entry in the oracle to an empty list
    updateLTUParameters([],[],[],[], fishIndex)

def resetMultifishOracle():
    # set parameters for each entry in the oracle to an empty list but NOT delete them
    for keys in multifishOracle.keys():
        resetRefFrameHistory(keys)

def clearMultifishOracle():
    # delete every entry in the oracle EXCPEPT the "0" entry. 
    # The parameters for this entry are just reset to empty lists.
    for keys in multifishOracle.keys():
        removeFishFromOracle(keys)
