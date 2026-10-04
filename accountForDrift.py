## Imports
import numpy as np

from pprint import pprint
import time
import math
from copy import copy

def matchFrames(seq1,seq2,drift):
    """Return overlapping crops of arrays shaped (frame, y, x).

    drift is (dx, dy): displacement of seq2 relative to seq1 in pixels.
    Positive dx/dy means the object moved right/down in seq2; cropping
    aligns it by selecting pixels farther right/down in seq2.
    Both returned arrays retain (frame, y, x) axis order.
    This XY drift convention differs from the legacy shifts.py YX convention.
    """
    dx = drift[0]
    dy = drift[1]

    # apply shifts
    # limit the drift that can be applied, to ensure we will always have at least one(!) pixel in each dimension
    rectF = [0,seq1[0].shape[0],0,seq1[0].shape[1]] # y1, y2, x1, x2
    rect = [0,seq2[0].shape[0],0,seq2[0].shape[1]] # y1, y2, x1, x2

    dx = np.minimum(dx, seq1[0].shape[1] - 1)
    dx = np.maximum(dx, -(seq1[0].shape[1] - 1))
    dy = np.minimum(dy, seq1[0].shape[0] - 1)
    dy = np.maximum(dy, -(seq1[0].shape[0] - 1))

    if dy<=0:
        rectF[0] = -dy # Note: will be a positive number
        rect[1] = rect[1]+dy
    else:
        rectF[1] = rectF[1]-dy
        rect[0] = dy
    if dx<=0:
        rectF[2] = -dx # Note: will be a positive number
        rect[3] = rect[3]+dx
    else:
        rectF[3] = rectF[3]-dx
        rect[2] = +dx

    seq1 = seq1[:,rectF[0]:rectF[1],rectF[2]:rectF[3]]
    seq2 = seq2[:,rect[0]:rect[1],rect[2]:rect[3]]

    return seq1,seq2
