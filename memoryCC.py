'''Modules for maintaining phase lock during time-lapse imaging.
Based on inter-frame cross-correlation between reference sets.
Each new reference frame sequence should be processed after determination.
This module can correct for known drift between sequences.'''

# Python Imports
import numpy as np
from pprint import pprint
import sys
import warnings
# Local Imports
import simpleCC as scc
import accountForDrift as afd
import shifts_global_solution as sgs
import plist_wrapper as jPlist
from image_loading import *
#from image_saving import *
from shifts import *
from shifts_global_solution import *

def processNewReferenceSequence(rawRefFrames,
                                thisPeriod,
                                thisDrift,
                                resampledSequences,
                                periodHistory,
                                driftHistory,
                                shifts,
                                knownPhaseIndex=0,
                                knownPhase=0,
                                numSamplesPerPeriod=80,
                                maxOffsetToConsider=2,
                                log=True):
    ''' Adapted from j_postacquisition.maintain_ref_frame_alignment

    Inputs:
    * rawRefFrames: a PxMxN numpy array representing the new reference frames
      (or a list of numpy arrays representing the new reference frames)
    * thisPeriod: the period for rawRefFrames (caller must determine the period)
    * thisDrift: the drift for rawRefFrames, in (x,y) order (caller must determine the drift).
      * if None no drift correction is used
    * resampledSequences: a list of resampled, previous reference frames
    * periodHistory: a list of the previous periods for resampledSequences
    * driftHistory: a list of the previous drifts for resampledSequences (in x,y order)
      * if no drift correction is used, this is a dummy variable
    * shifts: a list of shifts previously calculated for resampledSequences
    * knownPhaseIndex: the index of resampledSequences for which knownPhase applies
    * knownPhase: the phase (index) we are trying to match in knownPhaseIndex
    * numSamplesPerPeriod: the number of samples to use in resampling
    * maxOffsetToConsider: how far apart historically to make comparisons
      * should be used to prevent comparing sequences that are far apart and have little similarity

    Outputs:
    * resampledSequences: updated list of resampled reference frames
    * periodHistory: updated list of the periods for resampledSequences
    * driftHistory: updated list of the drifts for resampledSequences
      * if no drift correction is used, this is a dummy variable
    * shifts: updated list of shifts calculated for resampledSequences
    * globalShiftSolution[-1]: roll factor for latest reference frames
    * residuals: residuals on least squares solution'''

    # Deal with rawRefFrames type
    if type(rawRefFrames) is list:
        rawRefFrames = np.vstack(rawRefFrames)

    # Check that the reference frames have a consistent shape
    for f in range(1, len(rawRefFrames)):
        if rawRefFrames[0].shape != rawRefFrames[f].shape:
            # There is a shape mismatch.
            if log:
                # Return an error message and code to indicate the problem.
                print('Error: There is shape mismatch within the new reference frames. Frame 0: {0}; Frame {1}: {2}'.format(rawRefFrames[0].shape, f, rawRefFrames[f].shape))
            sys.stdout.flush()
            return (resampledSequences,
                    periodHistory,
                    driftHistory,
                    shifts,
                    -1000.0,
                    None)
    # And that shape is compatible with the history that we already have
    if len(resampledSequences) >= 1:
        if rawRefFrames[0].shape != resampledSequences[0][0].shape:
            # There is a shape mismatch.
            if log:
                # Return an error message and code to indicate the problem.
                print('Error: There is shape mismatch with historical reference frames. Old shape: {0}; New shape: {1}'.format(resampledSequences[0][0].shape, rawRefFrames[0].shape))
            sys.stdout.flush()
            return (resampledSequences,
                    periodHistory,
                    driftHistory,
                    shifts,
                    -1000.0,
                    None)

    if log:
        print('Using new reference frames with XY shape: {0}'.format(rawRefFrames[0].shape))

    # Add latest reference frames to our sequence set
    thisResampledSequence = scc.resampleImageSection(rawRefFrames,
                                                     thisPeriod,
                                                     numSamplesPerPeriod)
    resampledSequences.append(thisResampledSequence)
    periodHistory.append(thisPeriod)
    if thisDrift is not None:
        if len(driftHistory) > 0:
            # Accumulate the drift throughout history
            driftHistory.append([driftHistory[-1][0]+thisDrift[0],
                                 driftHistory[-1][1]+thisDrift[1]])
        else:
            driftHistory.append(thisDrift)
    else:
        if log:
            warnings.warn('No drift correction is being applied. This will seriously impact phase locking.',stacklevel=3)

    # Update our shifts array.
    # Compare the current sequence with recent previous ones
    if log:
        print(f'Drift history: {driftHistory}')
    if (len(resampledSequences) > 1):
        # Compare this new sequence against other recent ones
        firstOne = max(0, len(resampledSequences) - maxOffsetToConsider - 1)
        for i in range(firstOne, len(resampledSequences)-1):
            if log:
                print('---', i, len(resampledSequences)-1, '---')
            if thisDrift is None:
                alignment1, alignment2, rollFactor, score = scc.crossCorrelationRolling(resampledSequences[i][:numSamplesPerPeriod],
                                                                                        thisResampledSequence[:numSamplesPerPeriod],
                                                                                        numSamplesPerPeriod,
                                                                                        numSamplesPerPeriod)
            else:
                # apply drift correction first
                drift = [driftHistory[-1][0] - driftHistory[i][0],
                         driftHistory[-1][1] - driftHistory[i][1]]
                seq1, seq2 = afd.matchFrames(resampledSequences[i],
                                             resampledSequences[-1],
                                             drift)
                print(f'Drift correcting ({i}): {driftHistory[-1]}, {driftHistory[i]}. New shapes {seq1.shape}, {seq2.shape}')
                alignment1, alignment2, rollFactor, score = scc.crossCorrelationRolling(seq1,
                                                                                        seq2,
                                                                                        numSamplesPerPeriod,
                                                                                        numSamplesPerPeriod)
            shifts.append((i,
                           len(resampledSequences)-1,
                           rollFactor % numSamplesPerPeriod,
                           score))

    if log:
        pprint(shifts)

    # Linear regression for making historical shifts self consistent
    try:
        (globalShiftSolution, adjustedShifts, adjacentSolution, residuals, initialAdjacentResiduals) = sgs.MakeShiftsSelfConsistent(shifts,
                                                                                                                                    len(resampledSequences),
                                                                                                                                    numSamplesPerPeriod,
                                                                                                                                    knownPhaseIndex,
                                                                                                                                    knownPhase,
                                                                                                                                    log)
    except:
        print('Exception occurred during MakeShiftsSelfConsistent()')
        print(f'Input shifts {shifts}')
        print(f'Params {len(resampledSequences)}, {numSamplesPerPeriod}, {knownPhaseIndex}, {knownPhase}, {log}')
        sys.stdout.flush()      # Flush stdout so we can see all the information that led up to the point it went wrong
        raise
        
    if log:
        print('solution:')
        pprint(globalShiftSolution)

    # Calculate the residuals on the final solution - this is primarily useful for debugging
    residuals = np.zeros([len(globalShiftSolution), ])
    for i in range(len(globalShiftSolution)-1):
        for shift in shifts:
            if shift[1] == shifts[-1][1] and shift[0] == i:
                residuals[i] = globalShiftSolution[-1] - globalShiftSolution[i] - shift[2]
                break
        while residuals[i] > (numSamplesPerPeriod/2):
            residuals[i] = residuals[i]-numSamplesPerPeriod
        while residuals[i] < -(numSamplesPerPeriod/2):
            residuals[i] = residuals[i]+numSamplesPerPeriod

    if log:
        print('residuals:')
        pprint(residuals)
        print('Reference Frame rolling by: {0}'.format(globalShiftSolution[-1]))

    # Ensure correct sequencing between python print calls and printf calls from the C code that is calling us
    sys.stdout.flush()

    # Note for developers:
    # there are two other return statements in this function
    return (resampledSequences,
            periodHistory,
            driftHistory,
            shifts,
            globalShiftSolution[-1],
            residuals)

def RoIForReferenceHistory(resampledSequences):
    # This fuction has been relocated from j_postacquisition/maintain_ref_frame_alignment given this is the only function called by the LTU
    # Return the shape of the reference history.
    # We presume all have the same size (caller really should ensure this, or we will run into major problems!)
    if (len(resampledSequences) == 0):
        return (-1, -1)
    return resampledSequences[0][0].shape

def trimLTUHistory(resampledSequences,
                    periodHistory,
                    driftHistory,
                    shifts,
                    trimToLength):
    assert(len(resampledSequences) >= trimToLength)
    print(f"Trimming from initial sequence length {len(resampledSequences)} ({len(shifts)} shifts)")
    print(shifts)
    resampledSequences = resampledSequences[:trimToLength]
    periodHistory = periodHistory[:trimToLength]
    driftHistory = driftHistory[:trimToLength]
    shifts = shifts.copy()
    for n in range(len(shifts))[::-1]:
        i,j,_,_ = shifts[n]
        if j >= trimToLength:
            del shifts[n]
    print(f"Trimmed to sequence length {len(resampledSequences)} ({len(shifts)} shifts)")
    print(shifts)
    # Ensure correct sequencing between python print calls and printf calls from the C code that is calling us
    sys.stdout.flush()
    return resampledSequences,periodHistory,driftHistory,shifts

if __name__ == '__main__':
    print('Running toy example...This is STILL BROKEN')
    numStacks = 10
    stackLength = 10
    width = 10
    height = 10

    resampledSequencesDrift = []
    periodHistoryDrift = []
    shiftsDrift = []
    driftHistoryDrift = []

    resampledSequences = []
    periodHistory = []
    shifts = []
    driftHistory = []

    for i in range(numStacks):
        print('Stack {0}'.format(i))
        # Make new toy sequence
        thisPeriod = stackLength-0.5
        seq1 = (np.arange(stackLength)+np.random.randint(0, stackLength+1)) % thisPeriod
        print('New Sequence: {0}; Period: {1} ({2})'.format(seq1,
                                                            thisPeriod,
                                                            len(seq1)))
        seq2 = np.asarray(seq1, 'uint8').reshape([len(seq1), 1, 1])
        seq2 = np.repeat(np.repeat(seq2, width, 1), height, 2)

        # Run MCC without Drift
        resampledSequences, periodHistory, driftHistory, shifts, rollFactor, residuals = processNewReferenceSequence(seq2, thisPeriod, None, resampledSequences, periodHistory, driftHistory, shifts, knownPhaseIndex=0, knownPhase=0, numSamplesPerPeriod=80, maxOffsetToConsider=3, log=True)

        # Outputs for toy examples
        seqOut = (seq1+rollFactor) % thisPeriod
        print('Aligned Sequence (wout Drift): {0}'.format(seqOut))

        # Run MCC with Drift of [0,0]
        resampledSequencesDrift, periodHistoryDrift, driftHistoryDrift, shiftsDrift, rollFactor, residuals = processNewReferenceSequence(seq2, thisPeriod, [0, 0], resampledSequencesDrift, periodHistoryDrift, driftHistoryDrift, shiftsDrift, knownPhaseIndex=0, knownPhase=0, numSamplesPerPeriod=80, maxOffsetToConsider=3, log=True)

        # Outputs for toy examples
        seqOut = (seq1+rollFactor) % thisPeriod
        print('Aligned Sequence (with Drift): {0}'.format(seqOut))
        print('---')
