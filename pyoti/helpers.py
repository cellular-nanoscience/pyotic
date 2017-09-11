# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:49:29 2014

@author: Tobias Jachowski
"""
import collections
import io
import os
import numpy as np
try:
    import pandas as pd
    __pd__ = True
except ImportError:
    __pd__ = False
import sys
from scipy.ndimage import convolve1d
from scipy.ndimage import median_filter


if sys.platform == "linux" or sys.platform == "linux2":
    platform = 'linux'
elif sys.platform == "darwin":
    platform = 'darwin'
elif sys.platform == "win32":
    platform = 'win32'


def get_png_image(figure):
    figure.canvas.draw()
    buf = io.BytesIO()
    figure.savefig(buf, format='png')
    png = buf.getvalue()
    buf.close()
    return png


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


def moving_filter(data, window, moving_filter='mean', mode='reflect', cval=0.0,
                  origin=0):
    """
    Apply a moving filter to data.

    Parameters
    ----------
    data : numpy.ndarray
        The data to be filterd.
    window : int
        The window size of the moving filter.
    moving_filter : str, optional
        The filter to be used for the moving filter. Can be one of 'mean' or
        'median'.
    mode : str, optional
        mode       |   Ext   |         Data           |   Ext
        -----------+---------+------------------------+---------
        'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
        'reflect'  | 3  2  1 | 1  2  3  4  5  6  7  8 | 8  7  6
        'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
        'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
        'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3
        See 'scipy.ndimage.convolve1d' or 'scipy.ndimage.median_filter'
    cval : float, optional
        See 'scipy.ndimage.convolve1d' or 'scipy.ndimage.median_filter'
    origin : int, optional
        See 'scipy.ndimage.convolve1d' or 'scipy.ndimage.median_filter'
    """
    mode = mode or 'reflect'
    cval = cval or 0.0
    origin = origin or 0
    if moving_filter == 'mean' or moving_filter == 'average':
        return movingmean(data, window, mode=mode, cval=cval, origin=origin)
    else:  # if moving == 'median'
        return movingmedian(data, window, mode=mode, cval=cval, origin=origin)


def movingmean(data, window, mode='reflect', cval=0.0, origin=0):
    weights = np.repeat(1.0, window)/window
    # sma = np.zeros((data.shape[0] - window + 1, data.shape[1]))
    sma = convolve1d(data, weights, axis=0, mode=mode, cval=cval,
                     origin=origin)
    return sma


def movingmedian(data, window, mode='reflect', cval=0.0, origin=0):
    if data.ndim == 1:
        size = window
    else:
        size = (window, 1)
    smm = median_filter(data, size=size, mode=mode, cval=cval, origin=origin)
    return smm


def moving_mean(data, window):
    """
    Calculate a filtered signal by using a moving mean.

    This function is faster than the function `movingmean()`, but does not
    handle the lost edges of the data.

    Parameters
    ----------
    data : 1D numpy.ndarray of type float
        Data to calculate the rolling mean from.
    window : int
        Length of the window to calculate the rolling mean with.

    Returns
    -------
    1D numpy.ndarray of type float
        The data filtered with a rolling mean.
    """
    cumsum = np.cumsum(np.insert(data, 0, 0))
    return (cumsum[window:] - cumsum[:-window]) / window


def moving_mean_pandas(data, window):
    """
    Calculate a filtered signal by using a moving mean.

    Parameters
    ----------
    data : 1D numpy.ndarray of type float
        Data to calculate the rolling mean from.
    window : int
        Length of the window to calculate the rolling mean with.

    Returns
    -------
    1D numpy.ndarray of type float
        The data filtered with a rolling mean.
    """
    data = pd.Series(data)
    r = data.rolling(window=window)
    return r.mean()[window - 1:].get_values()


def calculate_means(data, samples=None, stds=False):
    """
    Calculate the means of `data`. If `samples_idx` is given, calculate the
    means of all samples contained in `samples_idx`.
    Additionally, the standard deviations can be returned, if stds is set to
    True.

    Parameters
    ----------
    data : numpy.ndarray
    samples : int, slice, or Iterable of int, slice, or index arrays
        Samples to calcualte the means from.
    stds : bool, optional
        Return also the standard deviations for the calculated means.

    Returns
    -------
    numpy.ndarray
        The calculated means
    (numpy.ndarray, numpy.ndarray)
        If `stds` is True, the calculated means and standard deviations.
    """
    if samples is None:
        samples = [slice(0, len(data))]
    if not isinstance(samples, collections.Iterable):
        samples = [samples]

    means = np.array([data[s].mean(axis=0) for s in samples], ndmin=2)

    if stds:
        stds = np.array([data[s].std(axis=0) for s in samples], ndmin=2)
        return means, stds

    return means


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


def files(directory, prefix=None, suffix=None, extension=None, sort_key=None):
    """
    Get filenames of a directory in the order sorted to their filename or a
    given key function.

    Parameters
    ----------
    directory : str
        The directory the files are located in.
    prefix : str
        Get only the files beginning with `prefix`.
    suffix : str
        Get only the files ending with `suffix`.
    extension : str, optional
        The extension of the files that should be returned. Default is
        '.txt'.
    sort_key : function
        Function to be applied to every filename found, before sorting.
    """
    prefix = prefix or ''
    suffix = suffix or ''
    extension = extension or ''
    files = [file_and_dir(filename=name, directory=directory)[2]
             for name in os.listdir(directory)
             if os.path.isfile(os.path.join(directory, name))
             and name.startswith(prefix)
             and name.endswith(''.join((suffix, extension)))]
    files.sort(key=sort_key)
    return files
