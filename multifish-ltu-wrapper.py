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

Each fish's inner dictionary contains the following keys.
Sequence indices are zero-based positions in resampledSequences.
All target phases and phase offsets are in resampled units (numSamplesPerPeriod samples per heartbeat),
rather than raw camera-frame units.

resampledSequences:
    Reference image sequence, representing a single heartbeat resampled with
    numSamplesPerPeriod samples. These are the sequences to be aligned.
periodHistory:
    The original heartbeat period, in camera frames, for each reference sequence.
    The alignment solver operates on the common resampled cycle length rather than these individual periods.
driftHistory:
    Accumulated XY drift in pixels, in (x, y) order, for each reference sequence.
    The difference between two entries is used to apply drift correction (if enabled)
    during the sequence-alignment process.
shifts:
    Pairwise alignment constraints (i, j, shift, score), normally with i < j.
    shift is the measured phase difference from sequence i to sequence j,
    wrapped to the resampled cycle. score is a metric of image mismatch; lower scores
    receive greater weight in the solver. Together these constraints are solved
    to determine the target phases across the history.
knownPhaseIndex, knownPhase:
    The reference sequence and target phase that anchor the alignment solution.
    If the user manually selects a new reference sequence this anchors the most recent sequence;
    the initial suggested target is also recorded this way.
    knownPhaseIndex starts at -1, which selects the latest sequence while the first reference is established.
lastStackStartIndex:
    The most recent reference successfully acquired by the standard refresh at
    the start of a z stack.
    When doing a z scan with a standard brightfield view (with changing focus),
    e.g. on a standard confocal microscope, we refresh the sync at regular intervals
    during a z scan. Then once the scan is complete we trim the history back to
    lastStackStartIndex
    None means no such boundary has yet been recorded, so a trim request would leave all history intact.
stackStartPhase:
    The equivalent target phase for lastStackStartIndex, cached from each refresh's
    existing solution and updated when the target is explicitly changed.
    We keep track of this for the unusual case where trimming removes the current anchor sequence
    (this should only happen if the user manually changes the target phase part-way through a z scan).
    If that happens, we apply this as the new knownPhase, to replace the one we have lost through trimming.
stackStartPhaseOffset:
    The solved phase at lastStackStartIndex minus the solved phase at the most
    recent sequence, modulo numSamplesPerPeriod. Adding it to a newly selected
    target on the most recent sequence updates stackStartPhase without another
    alignment solve being required. After trimming, both sequences are the boundary, so it is 0.
    Both cached phase values are None until a stack-start boundary is recorded.

The oracle retains only these two cached phase scalars, not the full solution
array. Clearing a fish's history clears its boundary, anchor and caches together.
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
              'knownPhase' : 0,
              'lastStackStartIndex' : None,
              'stackStartPhase' : None,
              'stackStartPhaseOffset' : None }

# The nested dict / oracle containing the LTU parameters for multiple fish.
# At startup we will not have any entries, but at least one should be added by the Spim GUI during its own startup
multifishOracle = dict()


def oracleDashboardRows():
    """Return a small, read-only snapshot for the helper's oracle dashboard.

    Keep this deliberately separate from get4LTUParameters/get6LTUParameters:
    those accessors create missing fish as a recovery measure, whereas merely
    displaying the dashboard must never change the oracle.

    The returned values are limited to Python ints and None so the Objective-C
    bridge never needs to inspect or copy any of the reference-frame arrays.
    A corrupt field affects only its own cell in the dashboard.
    """
    if type(multifishOracle) is not dict:
        raise TypeError('multifishOracle is not a dictionary')

    rows = []
    for uniqueFishID, parameters in list(dict.items(multifishOracle)):
        displayFishID = uniqueFishID if type(uniqueFishID) is int else None

        numReferenceSequences = None
        knownPhaseIndex = None
        if type(parameters) is dict:
            resampledSequences = dict.get(parameters, 'resampledSequences', [])
            if type(resampledSequences) in (list, tuple):
                numReferenceSequences = len(resampledSequences)

            candidateKnownPhaseIndex = dict.get(parameters, 'knownPhaseIndex', None)
            if type(candidateKnownPhaseIndex) is int:
                knownPhaseIndex = candidateKnownPhaseIndex

        rows.append((displayFishID, numReferenceSequences, knownPhaseIndex))

    # Fish IDs are normally integers. Any malformed-key rows remain visible at
    # the end of the table with an unavailable marker in the ID column.
    rows.sort(key=lambda row: (row[0] is None, row[0] if row[0] is not None else 0))
    return rows

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
        print(f'Updating knownPhase for unique fish ID {uniqueFishID} to fraction {fractionThroughSequence} (val {fractionThroughSequence * numSamplesPerPeriod})')
        parameters = multifishOracle[uniqueFishID]
        parameters['knownPhaseIndex'] = len(parameters['resampledSequences']) - 1
        parameters['knownPhase'] = fractionThroughSequence * numSamplesPerPeriod
        if parameters['lastStackStartIndex'] is not None:
            # This selection changes the target on the most recent reference,
            # which may later be discarded as an in-stack refresh. Preserve the
            # equivalent target in the retained stack-start reference as well.
            # The latest alignment already supplied the relative phase offset
            # between those sequences; changing the target does not change that
            # offset, so adding it gives the new fallback without another solve.
            # This updates only the cache: the selected sequence remains the
            # anchor unless a later trim actually removes it.
            parameters['stackStartPhase'] = (parameters['knownPhase'] + parameters['stackStartPhaseOffset']) % numSamplesPerPeriod
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

def getFractionalPhaseByAligningReferenceSequence(rawFrames, thisPeriod, thisDrift, maxOffsetToConsider, uniqueFishID, stackStart=False):
    # rawFrames is (frame, y, x); thisDrift is the (dx, dy) pixel displacement
    # from the preceding reference sequence (positive right/down).
    # Return target phase as a fraction of a cycle, or -1000.0 for shape mismatch.
    print(f'getFractionalPhaseByAligningReferenceSequence for unique fish ID {uniqueFishID}')
    ltuParameters = get6LTUParameters(uniqueFishID)
    resampledSequences, periodHistory, driftHistory, shifts, solution, _ = mcc.processNewReferenceSequence(rawFrames, thisPeriod, thisDrift, *ltuParameters, numSamplesPerPeriod, maxOffsetToConsider)
    if solution is None:   # shape mismatch: leave the boundary and caches intact
        sys.stdout.flush()
        return -1000.0
    shiftSolution = float(solution[-1])
    print(f'getFractionalPhaseByAligningReferenceSequence completed for unique fish ID {uniqueFishID} (result {shiftSolution:.3f}, frac {(shiftSolution/numSamplesPerPeriod)%1.0:.3f})')
    updateLTUParameters(resampledSequences, periodHistory, driftHistory, shifts, uniqueFishID)
    parameters = multifishOracle[uniqueFishID]
    if stackStart:
        # The acquisition notification identifies this particular refresh as a
        # stack start. Only record its new history index after alignment has
        # succeeded; ordinary/manual refreshes must not move the trim boundary.
        parameters['lastStackStartIndex'] = len(resampledSequences) - 1
    boundary = parameters['lastStackStartIndex']
    if boundary is not None:
        # Every successful refresh refines the equivalent target phase at the
        # boundary. Cache it now, while the full solution is available, for use
        # if a future trim discards the anchored sequence. Also cache the phase
        # difference to the newest sequence so a later explicit target selection
        # can update the fallback by scalar arithmetic. Neither operation changes
        # knownPhaseIndex/knownPhase, and no additional alignment solve is needed.
        # For a newly recorded boundary, boundary is the newest index and its
        # offset is zero. With no boundary, trimming is a no-op and no cache is needed.
        parameters['stackStartPhase'] = float(solution[boundary]) % numSamplesPerPeriod
        parameters['stackStartPhaseOffset'] = (float(solution[boundary]) - shiftSolution) % numSamplesPerPeriod
    # Note that we never actually use the residuals that get returned.
    # Only the newest solved phase is returned to the LTU helper app. The full
    # solution is temporary; the oracle retains only the boundary and two scalars.
    sys.stdout.flush()
    return (shiftSolution / numSamplesPerPeriod) % 1.0

def trimLTUHistory(uniqueFishID):
    print(f'Trim LTU history for fish {uniqueFishID}')
    ltuParameters = get4LTUParameters(uniqueFishID)
    parameters = multifishOracle[uniqueFishID]
    boundary = parameters['lastStackStartIndex']
    if boundary is None:
        print(f'Warning: no stack-start reference recorded for fish {uniqueFishID}; leaving LTU history unchanged', flush=True)
        return
    if not 0 <= boundary < len(ltuParameters[0]):
        raise ValueError(f'Invalid stack-start reference index {boundary} for fish {uniqueFishID}')
    trimToLength = boundary + 1
    anchorWillBeRemoved = parameters['knownPhaseIndex'] >= trimToLength
    # Normal calls cannot leave a recorded boundary without a cached phase:
    # the successful refresh that records it fills both caches in the same call,
    # and resetting history clears all three together. Guard against inconsistent
    # state from an incomplete/failed update before discarding the anchor.
    if anchorWillBeRemoved and parameters['stackStartPhase'] is None:
        raise ValueError(f'No cached stack-start target phase for fish {uniqueFishID}')
    returnTuple = mcc.trimLTUHistory(*ltuParameters, trimToLength)
    updateLTUParameters(*returnTuple, uniqueFishID)
    if anchorWillBeRemoved:
        print(f"Warning: trimming away phase anchor {parameters['knownPhaseIndex']} for fish {uniqueFishID}; transferring target to stack-start reference {boundary} at phase {parameters['stackStartPhase']}")
        parameters['knownPhaseIndex'] = boundary
        parameters['knownPhase'] = parameters['stackStartPhase']
    # Explicit target selections always refer to the most recent sequence in the
    # oracle. After trimming, that sequence is the stack-start boundary itself,
    # so a selection maps directly to stackStartPhase (no phase difference to add).
    # Discard the old offset to the removed latest sequence; the next successful
    # refresh will compute the offset to its newly appended sequence instead.
    parameters['stackStartPhaseOffset'] = 0.0
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
