'''Modules for phase matching two sequences based on cross-correlation.
This module includes all necessary functions.'''

# Python Imports
import numpy as np
from scipy.interpolate import interpn


def threePointTriangularMinimum(y1, y2, y3):
    # Fit an even V to three points at x=-1, x=0 and x=+1
    if y1 > y3:
        x = 0.5 * (y1-y3)/(y1-y2)
        y = y2 - x * (y1-y2)
    else:
        x = 0.5 * (y1-y3)/(y3-y2)
        y = y2 + x * (y3-y2)

    return x, y

def crossCorrelationScores(seq1,
                           seq2):
    '''Calculate temporal cross-correlation scores for (frame, flattened pixels) arrays.'''
    # Calculate cross-correlation from JT codes
    temp = np.conj(np.fft.fft(seq1, axis=0)) * np.fft.fft(seq2, axis=0)
    temp2 = np.fft.ifft(temp, axis=0)

    scores = np.sum(seq1*seq1) + np.sum(seq2*seq2) - 2 * np.sum(np.real(temp2),
                                                                axis=1)

    return scores


def minimumScores(scores,
                  useV=True):
    '''Calculates the minimum position and value in a list of scores.
    Can use V-fitting for sub-integer accuracy (useV=True).'''
    if useV:
        # V-fitting for sub-integer interpolation
        # Note that scores is a ring vector
        # i.e. scores[0] is adjacent to scores[-1]
        y1 = scores[np.argmin(scores)-1]  # Works even when minimum is at 0
        y2 = scores[np.argmin(scores)]
        y3 = scores[(np.argmin(scores)+1) % len(scores)]
        minPos, minVal = threePointTriangularMinimum(y1, y2, y3)
        minPos = (minPos + np.argmin(scores)) % len(scores)
    else:
        # Just use integer minimum
        minPos = np.argmin(scores)
        minVal = np.min(scores)

    return minPos, minVal


def matchSequenceSlicing(seq1,
                         seq2):
    '''Take two sequences and resample them to match the longer sequence.'''
    newLength = max(len(seq1), len(seq2))

    # Code review note (Oct 2026): the variable names x/y below denote array axes 1/2, actually image y/x!
    # Interpolation preserves the input (frame, y, x) layout without transposing.
    x = np.arange(0.0, seq1.shape[1])
    y = np.arange(0.0, seq1.shape[2])
    z1 = np.arange(0.0, seq1.shape[0])
    z2 = np.arange(0.0, seq2.shape[0])

    zOut = np.linspace(0.0, seq1.shape[0]-1, newLength)
    interpPoints = np.asarray(np.meshgrid(zOut, x, y, indexing='ij'))
    interpPoints = np.rollaxis(interpPoints, 0, 4)
    seq1 = interpn((z1, x, y), seq1, interpPoints)

    zOut = np.linspace(0, seq2.shape[0]-1, newLength)
    interpPoints = np.asarray(np.meshgrid(zOut, x, y, indexing='ij'))
    interpPoints = np.rollaxis(interpPoints, 0, 4)
    seq2 = interpn((z2, x, y), seq2, interpPoints)

    return seq1, seq2


def resampleImageSection(seq1,
                         thisPeriod,
                         newLength):
    '''Modified version of j_postacquisition.periods.ResampleImageSection'''
    # Input and return arrays are (frame, y, x); only the frame axis is resampled.

    result = np.zeros([newLength, seq1.shape[1], seq1.shape[2]], 'float')
    for i in range(newLength):
        desiredPos = (i / float(newLength)) * thisPeriod
        beforePos = int(desiredPos)
        afterPos = beforePos + 1
        remainder = (desiredPos - beforePos) / float(afterPos - beforePos)

        before = seq1[beforePos]
        after = seq1[int(afterPos % thisPeriod)]
        result[i] = before * (1 - remainder) + after * remainder

    return result


def crossCorrelationRolling(seq1,
                            seq2,
                            period1,
                            period2,
                            useV=True,
                            log=False,
                            numSamplesPerPeriod=80,
                            target=0):
    '''Phase matching two sequences based on cross-correlation.'''
    if log == 'toy':
        # Outputs for toy examples
        strout1 = []
        strout2 = []
        for i in seq1:
            strout1.append(i[0, 0])
        for i in seq2:
            strout2.append(i[0, 0])
        print('Original sequence #1:\t{0}'.format(strout1))
        print('Original sequence #2:\t{0}'.format(strout2))

    origLen1 = len(seq1)
    origLen2 = len(seq2)

    if period1 != numSamplesPerPeriod:
        seq1 = resampleImageSection(seq1, period1, numSamplesPerPeriod)
    if period2 != numSamplesPerPeriod:
        seq2 = resampleImageSection(seq2, period2, numSamplesPerPeriod)

    if log == 'toy':
        # Outputs for toy examples
        strout1 = []
        strout2 = []
        for i in seq1:
            strout1.append(i[0, 0])
        for i in seq2:
            strout2.append(i[0, 0])
        print('Resliced sequence #1:\t{0}'.format(strout1))
        print('Resliced sequence #2:\t{0}'.format(strout2))

    seq1 = makeArrayFromSequence(seq1[:numSamplesPerPeriod])
    seq2 = makeArrayFromSequence(seq2[:numSamplesPerPeriod])
    scores = crossCorrelationScores(seq1, seq2)
    print(f"crossCorrelationRolling: scores {scores}. Original lengths {origLen1} and {origLen2}")
    _, minValNoFit = minimumScores(scores, False)
    rollFactor, minVal = minimumScores(scores, useV)
    rollFactor = (rollFactor/len(seq1))*origLen2

    alignment1 = (np.arange(0, origLen1)-rollFactor) % origLen1
    alignment2 = np.arange(0, origLen2)

    if log is not False:  # will work for True and 'toy'
        print('Alignment 1:\t{0}'.format(alignment1))
        print('Alignment 2:\t{0}'.format(alignment2))
        print('Rolled by {0}'.format(rollFactor))

    if minValNoFit < 1e-5:
        # This really should only happen if we are replaying a dataset repeatedly.
        # In that scenario, we should get a minimum score value of 0.0.
        # However FP rounding errors may make it negative - which causes havoc.
        # In this scenario, limit it to a small positive number
        # (we don't want exactly zero either, as that would cause problems with our 1/score weightings)
        print(f'Suspicious negative value {minValNoFit}. Clipping this to small positive')
        minValNoFit = 1e-5
    if minVal < 0.02*minValNoFit:
        # This is probably only an issue for synthetic or bloodless datasets but I have seen the score go negative,
        # which is catastrophic for the weighted least-squares fitting used in the LTU,
        # since 1/sqrt(minVal) is computed as part of our weighting!
        print(f'Uncannily good fit {minVal}. Applying floor as a fraction of {minValNoFit} to make sure it does not get tiny or negative')
        minVal = 0.02*minValNoFit
    # convert to list for consistency with cjn-sequence-alignment.nCascadingNW
    return (list(alignment1),
            list(alignment2),
            (target+rollFactor) % origLen2,
            minVal)


def makeArrayFromSequence(seq):
    '''Utility function used as part of the FFT processing.'''
    result = np.zeros((len(seq), seq.shape[1]*seq.shape[2]))
    for i in range(len(seq)):
        result[i] = np.array(seq[i], copy=False).flatten()
    return result


if __name__ == '__main__':
    print('Running toy example with')
    str1 = [0, 1, 2, 3, 4, 3, 2, 2, 0]
    str2 = [4, 3, 2, 1, 0, 0, 1, 1, 1, 2, 3, 4]
    print('Sequence #1: {0}'.format(str1))
    print('Sequence #2: {0}'.format(str2))
    seq1 = np.asarray(str1, 'uint8').reshape([len(str1), 1, 1])
    seq2 = np.asarray(str2, 'uint8').reshape([len(str2), 1, 1])
    seq1 = np.repeat(np.repeat(seq1, 10, 1), 5, 2)
    seq2 = np.repeat(np.repeat(seq2, 10, 1), 5, 2)
    period1 = len(str1)-1
    period2 = len(str2)-1

    alignment1, alignment2, rollFactor, score = crossCorrelationRolling(seq1,
                                                                        seq2,
                                                                        period1,
                                                                        period2,
                                                                        log='toy')

    # Outputs for toy examples
    strout1 = []
    strout2 = []
    for i in alignment1:
        if i < 0:
            strout1.append(-1)
        else:
            strout1.append(str1[int(round(i % period1))])
    for i in alignment2:
        if i < 0:
            strout2.append(-1)
        else:
            strout2.append(str2[int(round(i % period2))])
    print('Aligned Sequence #1: {0}'.format(strout1))
    print('Aligned Sequence #2: {0}'.format(strout2))
