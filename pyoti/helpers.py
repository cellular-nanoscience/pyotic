# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:49:29 2014

@author: Tobias Jachowski
"""
import os
import numpy as np


def slicify(index, length=-1):
    """
    Takes an index as an instance of a list, a tuple, an np.ndarray, an int, or
    anything else. In the case of a list, a tuple, or a np.ndarray return a
    slice if possible, else return an np.ndarray. In the case of an int, return
    slice(int, int+1). In the case of anything else, return anything else.
    """
    # Convert list and tuple in an np.ndarray
    if isinstance(index, list) or isinstance(index, tuple):
        index = np.array(index)
    if isinstance(index, np.ndarray):
        # No elements in the array
        if len(index) == 0:
            return slice(0, 0, 1)
        # One element in the array
        if len(index) == 1:
            return slice(index[0], index[-1] + 1, 1)
        # More elements in the array
        if len(index) > 1:
            diff = np.diff(index)
            # Accept only strictly monotonically indices
            if diff[0] != 0 and np.all(diff[:-1] == diff[1:]):
                step = diff[0]
                start = index[0]
                stop = index[-1] + np.sign(step)
                # If the last element is at index 0 and the step is negative,
                # try to figure out the negative stop index. This is due
                # to the fact that a slice object behaves different for
                # numpy.ndarray than for the function `range()`.
                if stop < 0 and length < 0:
                    stop = step * len(index) - 1
                if stop < 0 and length >= 0:
                    stop = - length - 1
                return slice(start, stop, step)
    elif isinstance(index, int):
        return slice(index, index + 1, 1)

    return index


def listify(trace, length=-1):
    """
    Takes a trace as an instance of int, slice, or list and returns a list
    """
    if isinstance(trace, int) or isinstance(trace, str):
        return [trace]
    if isinstance(trace, slice):
        if trace.stop < 0 and trace.step < 0 and length < 0:
            length = abs(trace.stop) - 1
        elif trace.stop < 0 and trace.step > 0 and length < 0:
            length = abs(trace.stop) - 1
        elif trace.stop >= 0 and trace.step < 0 and length < 0:
            length = trace.start + trace.stop + 1
        elif trace.stop >= 0 and trace.step > 0 and length < 0:
            length = trace.stop
        else:  # length >= 0:
            length = length
        start, stop, step = trace.indices(length)
        return list(range(start, stop, step))
    if isinstance(trace, np.ndarray):
        return trace.tolist()
    return trace


def missing_elements(list1, list2):
    idx_existent = np.in1d(list1, list2)
    idx_missing = np.logical_not(idx_existent)
    missing_elements = listify(np.array(list1)[idx_missing])
    return missing_elements


def overlap_index(list1, list2):
    return np.nonzero(np.in1d(list1, list2))[0]


def skip_index(list1, list2):
    return np.nonzero(np.logical_not(np.in1d(list1, list2)))[0]


def file_and_dir(filename=None, directory=None):
    filename = filename or ""
    fdir = os.path.dirname(filename)
    ffile = os.path.basename(filename)

    ddir = directory or "."

    if (ffile == "" or ffile == "." or ffile == ".."):
        directory = os.path.join(ddir, filename, "")
        absdir = os.path.realpath(directory)
        return None, absdir, None

    directory = os.path.join(ddir, fdir, "")
    absdir = os.path.realpath(directory)
    absfile = os.path.join(absdir, ffile)

    return ffile, absdir, absfile
