# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:09:07 2016

@author: Tobias Jachowski

Segments are a 2D numpy.ndarray, with a start and stop index: ([start, stop])
"""
import numpy as np
import warnings
from scipy.ndimage import maximum_filter1d
from scipy.ndimage import minimum_filter1d


def get_excited_signal(signals):
    # Determine the excited signal out of two signals
    amplitude = signals.max(axis=0) - signals.min(axis=0)
    excited_signal = amplitude.argmax()
    return excited_signal


def get_excited_sections(signals, resolution, compare_time=0.01,
                         threshold=0.005, min_duration=0.2, validate=True):
    """
    Determine the indices (segments or range), where the first signal has a
    greater amplitude than the second signal.

    Parameters
    ----------
    signals : 2D array
        2D array of two signals to compare amplitude of.
    resolution : int
        The resolution of the `signals` in number of datapoints / s.
    compare_time : float
        Time in s for size of moving window used to compare the amplitude of
        the signals to each other.
    threshold : float
        Value, the amplitude of a signal is compared to, if only one signal is
        provided.
    min_duration : float
        Minimum duration in s a contiguous segment has to have to be finally
        accepted as one. Only used, if `validate` is True.
    validate : bool
        If True, excited segments (and non excited segments) need to be at
        least `min_duration` s.

    Returns
    -------
    tuple of segments and ranges :
    (
        (x_idx_segments, y_idx_segments),
        (x_range, y_range)
    )
    """
    window = int(np.round(compare_time * resolution))
    if window < 2:
        warnings.warn("Had to increase the compare_time for calculating "
                      "amplitudes within signal! Probably, increase the "
                      "resolution.")
        # compare_time = 2 / resolution
    window = max(2, window)

    greater = greater_amplitude(signals, window, threshold=threshold)

    min_length = 1
    if validate:
        min_length = max(1, int(np.round(min_duration * resolution)))

    x_excited = get_contiguous_segments(greater, min_length_high=min_length,
                                        min_length_low=min_length)
    if signals.ndim >= 2 and signals.shape[1] >= 2:
        y_excited = get_contiguous_segments(~greater,
                                            min_length_high=min_length,
                                            min_length_low=min_length)
    else:
        y_excited = np.empty((0, 2), dtype=int)

    x_excited_range = idx_segments_to_range(x_excited)
    y_excited_range = idx_segments_to_range(y_excited)

    return (x_excited, y_excited), (x_excited_range, y_excited_range)


def contains_greater_segment(x_segments, y_segments):
    """
    Determine, whether the first 2D segments array, contains a greater segment
    than the second 2D segments array.

    Parameters
    ----------
    x_segments : 2D segments array
        2D array with two columns. The first column contains the start and the
        second column the corresponding stop indices.
    y_segments : 2D segments array
        2D array with two columns. The first column contains the start and the
        second column the corresponding stop indices.

    Returns
    -------
    bool
        True, if x_segments contains a greater segment than y_segments.
        Otherwise, False.

    Examples
    --------
    positions = view.get_data(traces='positionXY', copy=False)
    x, y = sn.get_excited_sections(positions, view.samplingrate)[0]
    x_is_excited = sn.contains_greater_segment(x, y)
    """
    max_x = 0
    for x in x_segments:
        max_x = max(max_x, x[1] - x[0])
    max_y = 0
    for y in y_segments:
        max_y = max(max_y, y[1] - y[0])
    return max_x > max_y


def greatest_amplitude(signals, size, threshold=0.005):
    """
    Return 1D array with indices of the signal with the greatest amplitude.
    If signals.ndim == 1 or signals.shape[1] == 1 returns an array of boolean
    values which are True, if amplitude of signal is > threshold.

    Parameters
    ----------
    signals : 2D array
        2D array of a signal for which the amplitudes should be compared. If
        only one signal is given, the amplitude is compared to threshold.
        Otherwise, the amplitude is compared to the other amplitudes of the
        other signals.
        First dimension datapoints, second dimension different signals.
    size : int
        Size of the window the amplitudes should be calculated from and
        compared to each other. Needs to be >= 2.
    threshold : float
        Value, the amplitude of a signal is compared to, if only one signal is
        provided.

    Returns
    -------
    1D numpy.ndarray dtype=int or dtype=bool
        Index of where amplitude of first signal (first axis) is greater than
        the other signals or True where amplitude is greater than the
        threshold.
    """
    if size <= 1:
        raise ValueError("Size should be at least 2, to calculate amplitudes "
                         "within the signals!")
    amplitude = (maximum_filter1d(signals, size, axis=0)
                 - minimum_filter1d(signals, size, axis=0))
    if signals.ndim == 1 or signals.shape[1] == 1:
        greatest = amplitude > threshold
        greatest.shape = -1
    else:
        greatest = np.argmax(amplitude, axis=1)

    return greatest


def greater_amplitude(signals, size, threshold=0.005):
    """
    Determine, where a signal has a greater amplitude than a given threshold,
    or alternatively, if given, other signals.

    Parameters
    ----------
    signals : 2D array
        2D array of a signal for which the amplitudes should be compared. If
        only one signal is given, the amplitude is compared to `threshold`.
        Otherwise, the amplitude is compared to the other amplitudes of the
        other signals.
        First dimension datapoints, second dimension different signals.
    size : int
        Size of the window the amplitudes should be calculated from and
        compared to each other. Needs to be >= 2.
    threshold : float
        Value, the amplitude of a signal is compared to, if only one signal is
        provided.

    Returns
    -------
    1D numpy.ndarray dtype=bool
        True where amplitude of first signal (first axis) is greater than the
        other signals or the threshold.
    """
    greater = greatest_amplitude(signals, size, threshold=threshold)
    # Convert the index arrays into boolean arrays, True, where first signal
    # has greatest ampitude
    if signals.ndim >= 2 and signals.shape[1] >= 2:
        greater = (greater == 0)
    return greater


def get_contiguous_segments(plateaus, min_distance_center=1, min_length_high=1,
                            min_length_low=1, check_center_distance_first=True,
                            fuse=True, validate_truncated=True):
    # original source derived from
    # http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-bool_array-in-a-numpy-array
    # modified and extended by Tobias Jachowski
    """
    Find contiguous True values (plateau) of the boolean array `plateaus`.
    Return a 2D array where the first column is the start index of the region
    and the second column is the stop index (segments).

    Parameters
    ----------
    plateaus : np.ndarray of type bool
    min_distance_center : int
        Minimum distance, the center of two successive plateaus need to be
        apart from each other. If the distance is shorter, the stop of the
        leading and the start of the following plateau are neglected,
        effectively fusing the two plateaus.
    min_length_high : int
        Minimum number of contiguous True values a plateau has to have to get
        detected as one. If the plateau ist shorter, it is assumed to be no
        plateau and is therefore neglected. The minimum length of the plateaus
        usually is checked, after the distance of the center was checked.
        See parameter `check_center_distance_first`.
    min_length_low : int
        Minimum number of contiguous False values a valley has to have to get
        detected as one. If the valley is shorter, it is assumed to be part of
        a greater plateau, consisting of the preceding plateau, the short
        valley itself, and the succeeding plateau and is therefore neglected.
        The length of valleys is checked after the minimum length of the
        plateaus was checked.
    check_center_distance_first : bool
        Check `min_distance_center` before `min_length_high`.
    fuse : bool
        Fuse plateaus with a too short center distance, or delete the ones
        following the other ones too close.
    validate_truncated : bool, optional
        Validate the length of plateaus starting or stopping at the very first
        [0] or last [-1] index of the plateaus array, i.e. are possibly
        truncated plateaus. If set to True, even these plateaus need to have at
        least `min_length_high` datapoints. If set to False, even truncated -
        and therefore possibly too short - plateaus are detected.

    Returns
    -------
    2D numpy.ndarray of type int
        Segments of plateaus with start indices in the first and stop indices
        in the second column.
    """

    if plateaus.size == 0:
        return np.empty((0, 2), dtype=int)

    # Find the indices of changes in "plateaus"
    d = np.diff(plateaus)
    idx, = d.nonzero()

    # We need to start things after the change in "plateaus". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    lock_first = False
    if plateaus[0]:
        # If the start of plateaus is True prepend a 0
        idx = np.r_[0, idx]
        lock_first = not validate_truncated

    lock_last = False
    if plateaus[-1]:
        # If the end of plateaus is True, append the length of the array
        idx = np.r_[idx, plateaus.size]
        lock_last = not validate_truncated

    # No contiguous segments detected
    if idx.size < 2:
        return np.empty((0, 2), dtype=int)

    idx_start = idx[:-1:2]
    idx_stop = idx[1::2]

    def check_center_distance(idx_start, idx_stop):
        if min_distance_center == 1:
            return idx_start, idx_stop
        # Iteratively check distance of the center of one plateau (start, stop)
        # to the following one, and either fuse them or delete the next one, if
        # the distance is too small
        center = (idx_start + (idx_stop - idx_start) / 2).tolist()
        start = idx_start.tolist()
        stop = idx_stop.tolist()
        i = 0
        num_plateaus = len(center)
        while i < num_plateaus - 1:
            if center[i + 1] - center[i] < min_distance_center:
                if fuse:
                    # Correct the stop of the leading plateau to be the one
                    # of the following
                    stop[i] = stop[i + 1]
                    # Correct the center position of the now bigger leading
                    # plateau
                    center[i] = start[i] + (stop[i] - start[i]) / 2
                # Delete the following plateau and its corresponding center
                start.pop(i + 1)
                stop.pop(i + 1)
                center.pop(i + 1)

                # Correct number for deleted plateau and go back one step to
                # check also the distance of the current plateau to the plateau
                # after the next one (now the next one)
                num_plateaus -= 1
                i -= 1
            i += 1
        return np.array(start), np.array(stop)

    def check_length_high(idx_start, idx_stop):
        if min_length_high == 1:
            return idx_start, idx_stop
        # Check length of plateaus first
        # (stop_high - start_high) or (uneven_index - even_index) >= min_length
        plateau = idx_stop - idx_start >= min_length_high
        # If stops are locked, make sure to reset them to True
        plateau[0] = plateau[0] or lock_first
        plateau[-1] = plateau[-1] or lock_last
        # Select only accepted starts/stops (plateaus)
        idx_start = idx_start[plateau]
        idx_stop = idx_stop[plateau]
        return idx_start, idx_stop

    checks = [check_length_high, check_center_distance]
    if check_center_distance_first:
        checks.reverse()
    for check in checks:
        idx_start, idx_stop = check(idx_start, idx_stop)

    # Only one plateau detected, ignore all valleys
    if idx_start.size <= 1:
        idx = np.sort(np.r_[idx_start, idx_stop])
        idx.shape = (-1, 2)
        return idx

    # Check length of valleys
    # (stop_low - start_low) or (even_index - uneven_index) >= min_length
    valley = np.where(idx_start[1:] - idx_stop[:-1] >= min_length_low)[0]
    # Select only accepted stops/starts (valleys) and keep first start and last
    # stop of plateaus
    start = idx_start[np.r_[0, valley + 1]]
    stop = idx_stop[np.r_[valley, idx_stop.size - 1]]

    # Concatenate and sort passed start indices
    idx = np.sort(np.r_[start, stop])
    # Reshape the result into two columns
    idx.shape = (-1, 2)

    return idx


def get_steps(signal, resolution, step_time=0.1, min_dwell_time=0.0):
    """
    Detect positive and negative steps in a signal.

    Parameters
    ----------
    signal : 1D numpy.ndarray
        Signal that contains steps to be detected.
    resolution : float
        Resolution of the signal in Hz.
    step_time : float
        The time a step takes in s. Should at least be greater as the greatest
        period that can be expected due to noise.
    min_dwell_time : float
        Min time between two succesive steps in s. If dwell time between two
        pre_detected steps is shorter the second step will be ignored.

    Returns
    -------
    tuple of 1D int index arrays
        (pos_steps, neg_steps)
    """
    window = max(1, int(step_time * resolution))
    min_dwell_length = max(1, int(min_dwell_time * resolution))
    return _get_steps(signal, window, min_dwell_length)


def _get_steps(signal, window, min_dwell_length=1):
    X = signal
    Xl = X[:-window]
    Xr = X[window:]

    pos_step_plateaus = Xl < Xr
    neg_step_plateaus = Xl > Xr

    # The detected plateaus caused by the steps should have a length of at
    # least size window (with a factor of safety due to noise)
    min_length_high = window * 0.95
    # The detected valleys caused by the dwell time (dwell length) are
    # shortened by a value of size window.
    min_length_low = max(1, (min_dwell_length - window) * 0.95)

    pos_steps = get_contiguous_segments(pos_step_plateaus,
                                        min_length_high=min_length_high,
                                        min_length_low=min_length_low,
                                        validate_truncated=True)
    neg_steps = get_contiguous_segments(neg_step_plateaus,
                                        min_length_high=min_length_high,
                                        min_length_low=min_length_low,
                                        validate_truncated=True)

    # The stop indices is the start of every new dwell plateau (directly
    # after the step)
    # pos_steps = pos_steps[:,1]
    # neg_steps = neg_steps[:,1]

    # The mean of start and stop is window/2 left of the middle of the step
    pos_steps = np.round(pos_steps.mean(axis=1) + window/2).astype(int)
    neg_steps = np.round(neg_steps.mean(axis=1) + window/2).astype(int)

    return pos_steps, neg_steps


def get_extrema(signal, resolution=1, highest_frequency=1,
                reduce_false_positives=False, idx_plateaus=False):
    """
    Determine the indices of Minima und Maxima of a signal.

    Parameters
    ----------
    signal : 1D numpy.ndarray
        The signal to be analysed.
    resolution : float
        Resolution of the signal (e.g. in Hz).
    highest_frequency : float
        Highest expected frequency in Hz of occurence of minima/maxima in the
        signals. If you set a too low frequency, extrema will not be detected.
        A frequency to high will result in to many falsly positive detected
        extrema. See also parameter `reduce_false_extrema`.
    reduce_false_positives : bool
        Set to True, to reduce falsly positive detected extrema. As a side
        effect, the extrema within half a wavelength of the given highest
        frequency at the ends of the signal could be shifted or even not be
        detected. Make sure to properly set `highest_frequency`, first.
        Another option to reduce False positives is to smoothen the signal,
        before using this function.
    idx_plateaus : bool
        If idx_plateaus is True, return not only the indices of the extrema,
        but also the start and the stop indices (segments) of the plateaus
        corresponding to the detected extrema.

    Returns
    -------
    tuple of numpy.ndarray of int
        Indices of detected minima and maxima (minima, maxima).
    """
    # The time should be adjusted according to the expected frequency and
    # should be < 1/2 of the period of the highest frequency expected (one half
    # wave equals to one extremum, factor=0.5). This time is used to calculate
    # a window size, for which the values of two datapoints are compared to be
    # greater or smaller -> extrema detection.
    # For reliable extremum detection at the ends of the signal, time should
    # even be less (factor=0.05).
    # Time <- depending on noise and amplitude, too!
    # -> higher frequency/amplitude, smaller time
    # greater time -> less false positives
    if reduce_false_positives:
        factor = 0.45  # < 50% of the period (95% of 50%)
    else:
        factor = 0.05  # <= 5% of the period
    time = factor/highest_frequency
    window = max(1, int(time * resolution))
    return _get_extrema(signal, window, idx_plateaus=idx_plateaus)


def _get_extrema(signal, window, idx_plateaus=False):
    """
    Get attributes of a signal like pattern X:

    Parameters
    ----------
    signal
    window : int
        should be greater than 50 points for A = 6 µm and f = 0.2 Hz, because
        the noise of the piezo feed back is such that it takes roughly 50
        datapoints to clearly discriminate between a data point with a true
        local maximum value and a data point with a local maximum due to noise
        (plateau of ~50) Should be smaller than datapoints for one half period,
        too. 100 is a good value for triangle signal 400 is a save value for
        sinusiodal signal 350 works for sinusiodal signal
    idx_plateaus : bool
        If idx_plateaus is True, return not only the indices of the extrema,
        but also the start and the stop indices (segments) of the plateaus
        corresponding to the detected extrema.

    Returns
    -------
    tuple : tuple of numpy.ndarrays of type int
        Tuple of arrays of minima and maxima.
    """
    # For future implementations of _get_extrema() or another implementation,
    # probably consider one of the solutions presented in
    # https://blog.ytotech.com/2015/11/01/findpeaks-in-python/
    X = signal

    Xl = X[:-2*window]
    Xm = X[window:-window]
    Xr = X[2*window:]

    # find indices of data points, which are smaller/greater than the
    # left/right ones:
    # left datapoint = comparepoint - window
    # right datapoint = comparepoint + window
    # This will create a plateau:
    minima_plateaus = np.logical_and(Xm < Xr, Xm < Xl)
    maxima_plateaus = np.logical_and(Xm > Xr, Xm > Xl)

    # Get start/stop indices of plateaus, make sure they have at least a length
    # of window and are at least a length of window apart, accept truncated
    # extrema plateaus, if they start or end at the very first [0] or last [-1]
    # index of the plateaus array.
    # Consider noise to shorten the length of plateaus and valleys to up to 5 %
    min_length_high = window * 0.95
    min_length_low = window * 0.95

    minima = get_contiguous_segments(minima_plateaus,
                                     min_length_high=min_length_high,
                                     min_length_low=min_length_low,
                                     validate_truncated=False)
    maxima = get_contiguous_segments(maxima_plateaus,
                                     min_length_high=min_length_high,
                                     min_length_low=min_length_low,
                                     validate_truncated=False)

    # correct the indices for the shift due to truncating Xm (Xl and Xr)
    minima += window
    maxima += window

    # calculate center position of start and stop indices
    minima_c = np.round(minima.mean(axis=1)).astype(int)
    maxima_c = np.round(maxima.mean(axis=1)).astype(int)

    if idx_plateaus:
        return minima_c, maxima_c, minima, maxima

    return minima_c, maxima_c


def get_value_segments(signal, value=0.0, std_noise=0.01, max_length=17,
                       min_length=1):
    """
    Search for the locations where the signal is equal to value and return
    index segments of these locations.

    Parameters
    ----------
    value : float, optional
        The value that should be looked up in `signal`.
    std_noise : float, optional
        Root mean square of the noise of the signal.
    max_length : int, optional
        Max number of contiguous datapoints that are allowed to be outside the
        `value` +/- `std_noise` range, but still be assumed to be equal to the
        `value`. On average, one out of 1 / (0.32^`max_length`) datapoints will
        falsly be detected as not beeing equal to the `value`. Assuming the
        `std_noise` was set correctly, this would be 1 out of ~ 258,000,000
        datapoints for the default value 17.
    min_length : int, optional
        Min number of contiguous datapoints that need to be within the range
        `value` +/- `std_noise`, to be detected as equal to `value`.

    """
    # Get the boolean idx array of the datapoints beeing within 0 +/-
    # std_noise.
    s = signal - value
    base_g = s > - std_noise
    base_s = s < std_noise
    base = np.logical_and(base_g, base_s)
    base = base.flatten()
    seg = get_contiguous_segments(base, min_length_low=max_length,
                                  min_length_high=min_length)
    return seg


def get_segment_centers(segments):
    if len(segments) == 0:
        return np.empty((0, 2))
    diff = segments[:, 1] - segments[:, 0]
    centers = np.round((segments[:, 0] + diff / 2)).astype(int)
    return centers


def center_sampled_segments(segments, centers=None, decimate=1, as_list=False,
                            slices=False):
    if as_list:
        idx = []
    else:
        if slices:
            dtype = object
        else:
            dtype = int
        idx = np.empty(0, dtype=dtype)

    if slices:
        f_idx = slice
    else:
        f_idx = np.arange

    if centers is None:
        centers = get_segment_centers(segments)

    idx = [f_idx(segment[0] + (center - segment[0]) % decimate or decimate,
                 segment[1],
                 decimate)
           for segment, center in zip(segments, centers)]

    if not as_list and not slices:
        idx = np.hstack(idx)

    return idx


def baseline_idx(signals, value=0.0, std_noise=0.01, max_length=17,
                 decimate=1):
    """
    Calculate the indices of where the `signals` are equal to `value` within
    a given standard deviation `std_noise`.
    This method is based on the function `signal.get_value_segments()` and
    `signal.center_sampled_segments()`.

    Parameters
    ----------
    signals : 2D numpy.ndarray of type float
        Array containing the signals in the columns.
    value : float, optional
        See function `pyoti.evaluate.signal.get_value_segments()`
    std_noise : float, optional
        See function `pyoti.evaluate.signal.get_value_segments()`
    max_length : int, optional
        See function `pyoti.evaluate.signal.get_value_segments()`
    decimate : int, optional
        See function `pyoti.evaluate.signal.center_sampled_segments()`
    """
    if signals.ndim == 1:
        signals = signals[:, np.newaxis]
    idx = np.empty(0, dtype=int)
    for signal in signals.T:
        segments = get_value_segments(signal, value=value, std_noise=std_noise,
                                      max_length=max_length)
        centers = get_segment_centers(segments)
        i = center_sampled_segments(segments, centers=centers,
                                    decimate=decimate)
        idx = np.r_[idx, i]
    idx.sort()
    return idx


def get_sections(signal, minima, maxima):
    """
    Separate a signal into rising (rise) and falling (fall) sections.
    And in turn separate these sections into positive (pos) and
    negative (neg) sections.

    Returns
    -------
    tuple of segments and ranges :
    (
        (posrise, posfall, negfall, negrise, pos, neg, rise, fall),
        (posrise_range, posfall_range, negfall_range, negrise_range,
         pos_range, neg_range, rise_range, fall_range)
    )
    """
    X = signal
    pos = X > 0
    neg = ~ pos

    # Get the indices of the change of sign.
    # The indices should have a difference of period/2 (except the difference
    # from first and to last datapoint).
    # A minimum length of period/3 should be save to exclude small local
    # sign changes, but still detect true (global) sign changes
    # Include the very first datapoint of pos/neg range and shift index by 1
    # min_period = get_period(minima, maxima, mode='min')
    # length = int(np.round(min_period/3))
    pos = get_contiguous_segments(pos)
    neg = get_contiguous_segments(neg)

    # Make sure the ranges are contiguous and not interrupted like the original
    # ones received by comparison of sign of X.
    # pos_range = np.hstack([ np.arange(change_of_sign[i],
    #                                   change_of_sign[i + 1])
    #                         for i in np.arange(0,
    #                                            change_of_sign.size - 1,
    #                                            2)])
    # neg_range = np.hstack([ np.arange(change_of_sign[i],
    #                                   change_of_sign[i + 1])
    #                         for i in np.arange(1,
    #                                            change_of_sign.size - 1,
    #                                            2)])
    pos_range = idx_segments_to_range(pos)
    neg_range = idx_segments_to_range(neg)

    # calculate all extrema, extend first and last extremum and avoid invalid
    # indices
    extrema = np.sort(np.r_[minima, maxima])

    if extrema.size == 0:
        rise_range = np.empty((0, 2), dtype=int)
        fall_range = np.empty((0, 2), dtype=int)
        rise = np.empty((0, 2), dtype=int)
        fall = np.empty((0, 2), dtype=int)
    else:
        if extrema.size == 1:
            # Only one maximum or minimum -> no info about period.
            # Take all indices before and after extremum
            first = 0
            last = X.size
        else:
            max_period = get_period(minima, maxima, mode='max')
            half_period = int(np.round(max_period / 2))
            first = max(0, extrema[0] - half_period)
            last = min(X.size, extrema[-1] + half_period)

        extrema = np.r_[first, extrema, last]

        # Assume very fist extremum is (appended) minimum:
        # maxima[0] < minima[0]
        # -> extremum[0] -> extremum[1] = rise
        if extrema.size % 2 == 0:
            idx_rise = slice(0, None)
            idx_fall = slice(1, -1)
        else:
            idx_rise = slice(0, -1)
            idx_fall = slice(1, None)

        rise = extrema[idx_rise].reshape(-1, 2)
        fall = extrema[idx_fall].reshape(-1, 2)

        # Very first extremum is (appended) maximum
        # exchange fall and rise
        if ((maxima.size == 1 and extrema.size == 1) or
           (maxima.size >= 1 and minima.size >= 1 and minima[0] < maxima[0])):
            rise, fall = fall, rise

        rise_range = idx_segments_to_range(rise)
        fall_range = idx_segments_to_range(fall)

    posrise_range = np.intersect1d(rise_range, pos_range)
    posfall_range = np.intersect1d(fall_range, pos_range)
    negfall_range = np.intersect1d(fall_range, neg_range)
    negrise_range = np.intersect1d(rise_range, neg_range)

    posrise = idx_range_to_segments(posrise_range)
    posfall = idx_range_to_segments(posfall_range)
    negfall = idx_range_to_segments(negfall_range)
    negrise = idx_range_to_segments(negrise_range)

    # posrise = []
    # for stop in np.intersect1d(rise, pos_range + 1):
    #     idx = int(np.floor(np.where(np.in1d(rise, stop))[0] / 2))
    #     start = pos_range[
    #                 np.min(np.intersect1d(np.where(pos_range <= stop)[0],
    #                              np.where(pos_range >= rise[idx, 0])[0]))]
    #     #if start != stop:
    #     posrise.append(np.array([start, stop]))
    # posrise = np.array(posrise)
    # posfall = []
    # for start in np.intersect1d(fall, pos_range):
    #     idx = int(np.floor(np.where(np.in1d(fall, start))[0] / 2))
    #     stop = pos_range[
    #                np.max(np.intersect1d(np.where(pos_range >= start)[0],
    #                             np.where(pos_range <= fall[idx, 1])[0]))] + 1
    #     #if start != stop:
    #     posfall.append(np.array([start, stop]))
    # posfall = np.array(posfall)
    # negfall = []
    # for stop in np.intersect1d(fall, neg_range + 1):
    #     idx = int(np.floor(np.where(np.in1d(fall, stop))[0] / 2))
    #     start = neg_range[
    #                 np.min(np.intersect1d(np.where(neg_range <= stop)[0],
    #                              np.where(neg_range >= fall[idx, 0])[0]))]
    #     #if start != stop:
    #     negfall.append(np.array([start, stop]))
    # negfall = np.array(negfall)
    # negrise = []
    # for start in np.intersect1d(rise, neg_range):
    #     idx = int(np.floor(np.where(np.in1d(rise, start))[0] / 2))
    #     stop = neg_range[
    #                np.max(np.intersect1d(np.where(neg_range >= start)[0],
    #                             np.where(neg_range <= rise[idx, 1])[0]))] + 1
    #     #if start != stop:
    #     negrise.append(np.array([start, stop]))
    # negrise = np.array(negrise)

    return ((posrise, posfall, negfall, negrise, pos, neg, rise, fall),
            (posrise_range, posfall_range, negfall_range, negrise_range,
             pos_range, neg_range, rise_range, fall_range)
            )


def get_period(minima, maxima, mode='mean'):
    """
    Calculate the period of given indices of extrema (minima and maxima) of an
    oscillation.

    Parameters
    ----------
    minima : 1D array
        1D array of indices of minima of an oscillation.
    maxima : 1D array
        1D array of indices of minima of an oscillation.
    mode : str
        Return the 'mean', 'max', 'min', or 'median' period.

    Returns
    -------
    float or int
        length of the period of the oscillation
    """
    extrema = np.sort(np.r_[minima, maxima])
    periods = np.diff(extrema) * 2

    if periods.size >= 1:  # at least 2 extrema
        if mode == 'mean' or mode == 'average':
            period = np.mean(periods)
        elif mode == 'max':
            period = np.max(periods)
        elif mode == 'min':
            period = np.min(periods)
        else:  # 'median'
            period = np.median(periods)
    else:
        period = 0

    return period


def idx_to_idx_segments(idx, start=None, stop=None):
    """
    Create stop indices
    Concatenate and sort start (and stop) indices
    """
    if idx.size >= 2:
        idx_segs = np.sort(np.r_[idx, idx[1:-1]])
    else:  # size is 1 or empty
        idx_segs = np.array([], dtype=int)

    if idx.size >= 2:
        if start is not None and start < idx[0]:
            idx_segs = np.r_[start, idx[0], idx_segs]
        if stop is not None and stop > idx[-1]:
            idx_segs = np.r_[idx_segs, idx[-1], stop]
    if idx.size == 1:
        if start is not None and start < idx[0]:
            idx_segs = np.r_[start, idx[0], idx_segs]
        if stop is not None and stop > idx[-1]:
            idx_segs = np.r_[idx_segs, idx[-1], stop]

    # Reshape the result into two columns
    idx_segs.shape = (-1, 2)
    return idx_segs


def idx_segments_to_range(segments, decimate=None):
    """
    Convert a 2D segments array into a 1D range array.
    """
    idx = np.empty(0, dtype=int)
    if segments.size > 0:
        idx = np.hstack([np.arange(segment[0], segment[1], decimate)
                        for segment in segments])
    return idx


def idx_segments_to_slices(segments, decimate=None):
    """
    Convert a 2D segments array into a 1D array of slices.
    """
    if decimate is not None:
        decimate = int(np.round(decimate))
    slices = []
    if segments.size > 0:
        slices = np.array([slice(segment[0], segment[1], decimate)
                          for segment in segments])
    return slices


def idx_range_to_segments(range):
    """
    Convert a 1D range array into a 2D segments array.
    """
    segments = np.empty((0, 2), dtype=int)
    if range.size >= 1:
        idx = np.where(np.diff(range) > 1)[0]
        start = np.r_[range[0], range[idx + 1]]
        stop = np.r_[range[idx], range[-1]] + 1
        segments = np.sort(np.r_[start, stop])
        segments.shape = (-1, 2)
    return segments


def get_min_size_segments(segments, size):
    segments = segments[segments[:, 1] - segments[:, 0] >= size]
    return segments


def shift_and_size_segments(segments, shift, size, min_start=-1, max_stop=0):
    segments = shift_segments(segments, shift)
    segments = size_segments(segments, size)

    segments = limit_segments(segments, min_start=min_start, max_stop=max_stop)

    return segments


def size_segments(segments, size, min_start=-1, max_stop=0):
    if size > 0:
        segments = np.sort(np.r_[segments[:, 0], segments[:, 0] + size])
        segments.shape = (-1, 2)

    segments = limit_segments(segments, min_start=min_start, max_stop=max_stop)

    return segments


def shift_segments(segments, shift, min_start=-1, max_stop=0):
    # shift ("rotate") the indices of segments
    segments = segments + shift

    segments = limit_segments(segments, min_start=min_start, max_stop=max_stop)

    return segments


def limit_segments(segments, min_start=-1, max_stop=0):
    # create copy to protect outer array
    segments = segments.copy()

    if segments.ndim == 1:
        # if no segments but indices
        if min_start > -1:
            segments[segments < min_start] = min_start
        if max_stop > 0:
            segments[segments > max_stop - 1] = max_stop - 1
    else:
        if min_start > -1:
            starts = segments[:, 0]
            stops = segments[:, 1]
            starts[starts < min_start] = min_start
            stops[stops < min_start + 1] = min_start + 1
        if max_stop > 0:
            starts = segments[:, 0]
            stops = segments[:, 1]
            starts[starts > max_stop - 1] = max_stop - 1
            stops[stops > max_stop] = max_stop

    return segments


def get_excited_offset(signal, minima, maxima, period):
    X = signal

    # Calculate voltage (offset) of signal
    offset_min_pos = minima + np.round(period/4)
    left_first_min_offset = minima[0] - np.round(period/4)
    if left_first_min_offset > 0:
        offset_min_pos = np.append(left_first_min_offset, offset_min_pos)
    offset_min_pos = offset_min_pos[offset_min_pos < X.shape[0]].astype(int)

    offset_max_pos = maxima + np.round(period/4)
    left_first_max_offset = maxima[0] - np.round(period/4)
    if left_first_max_offset > 0:
        offset_max_pos = np.append(left_first_max_offset, offset_max_pos)
    offset_max_pos = offset_max_pos[offset_max_pos < X.shape[0]].astype(int)

    offset = np.mean(np.append(X[offset_min_pos], X[offset_max_pos]))

    return offset


def outliers(signal, iqr_factor=1.5):
    # Filter out outliers
    # Points with mean intensities that are further than 1.5 times the IQR
    # away from 1st or 3rd quartile
    # iqr_factor = iqr_factor # standard 1.5, strong 3.0
    qrt1 = np.nanpercentile(signal, 25, interpolation='midpoint')
    qrt3 = np.nanpercentile(signal, 75, interpolation='midpoint')
    iqr = qrt3 - qrt1
    upper = qrt3 + iqr * iqr_factor
    lower = qrt1 - iqr * iqr_factor
    wo_outls = np.logical_and(signal <= upper, signal >= lower)
    w_outls = np.logical_or(signal > upper, signal < lower)
    return wo_outls, w_outls


def validate_extrema(extrema):
    # Rewrite this faulty function!!!
    return
    # ### exclude outliers
    # exclude outliers by calculating the quartiles
    # if the signal has a constant frequency the distance of one extremum to
    # the next one should be constant
    diff = extrema[1:] - extrema[:-1]
    median = np.median(diff)
    lower = 0.95 * median
    upper = 1.05 * median
    idx, = np.where(np.logical_and(diff > lower, diff < upper))
    if idx.size == 0:
        idx = slice(0, 1, 1)
    else:
        idx = np.append(idx, idx[-1] + 1)
    return extrema[idx]
