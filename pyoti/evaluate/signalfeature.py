# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 22:03:33 2016

@author: Tobias Jachowski
"""
import numpy as np
import warnings

from . import signal as sn


class SignalFeature(object):
    """
    A `SignalFeature` analyzes a signal and detects specific sections in the
    signal (like steps, plateaus, etc.). It offers properties to easily get the
    indices/segments of the different sections. It provides the method
    `update()` to update the detected sections. `SignalFeature` is meant to be
    subclassed. The subclassed classes have to provide the function
    `_update()`, which will do the actual update of the sections.
    """
    def __init__(self, signals, resolution, **kwargs):
        self.updated = False
        self.update(signals, resolution)

    def update(self, signals, resolution):
        self._update(signals, resolution)
        self.updated = True

    def _update(self, signals, resolution):
        # do the update here
        pass


class StepFinder(SignalFeature):
    """
    A `StepFinder` finds the indices of positive and negative steps in a
    signal.
    """
    def __init__(self, signals, resolution, step_time=0.1, min_dwell_time=0.0,
                 include_first_last=False, **kwargs):
        """
        Parameters
        ----------
        signals : 2D numpy.ndarray
            Signal that contains steps to be detected.
        resolution : float
            Resolution of the signals in Hz.
        step_time : float
            The time a step takes in s. Should at least be greater as the
            greatest period that can be expected due to noise.
        min_dwell_time : float
            Min time between two successive steps in s. If dwell time between
            two detected steps is shorter the second step will be ignored.
        include_first_last : bool
            If `include_first_last` is True, the attributes `plateau_segments`,
            `pos_plateau_segments` and `neg_plateau_segments` will include the
            first datapoint to the first step and the last step to the last
            datapoint as a plateau.

        Attributes
        ----------
        include_first_last
        step_time
        min_dwell_time
        steps : 1D numpy.ndarray dtype=int
        pos_steps : 1D numpy.ndarray dtype=int
        neg_steps : 1D numpy.ndarray dtype=int
        plateau_segments : 2D numpy.ndarray dtype=int
            The first dimension contains the plateaus, the second dimension the
            start and stop value.
        pos_plateau_segments : 2D numpy.ndarray dtype=int
            The first dimension contains the plateaus, the second dimension the
            start and stop value.
        neg_plaeau_segments : 2D numpy.ndarray dtype=int
            The first dimension contains the plateaus, the second dimension the
            start and stop value.
        """
        self._step_time = step_time
        self._min_dwell_time = min_dwell_time
        self._include_first_last = include_first_last

        super().__init__(signals=signals, resolution=resolution)

    def _update(self, signals, resolution):
        signal = signals[:, 0]
        self.max_stop = signal.size
        self._steps = {}  # 'pos' and 'neg' steps of signal
        self._steps['pos'], self._steps['neg'] = sn.get_steps(
            signal, resolution, self.step_time, self.min_dwell_time)

    @property
    def include_first_last(self):
        return self._include_first_last

    @include_first_last.setter
    def include_first_last(self, include_first_last):
        self._include_first_last = include_first_last
        self.updated = False

    @property
    def step_time(self):
        """
        step_time : (s) time a step is supposed to take
            should at least be greater as the greatest period to be expected
            due to noise
        """
        return self._step_time

    @step_time.setter
    def step_time(self, step_time):
        self._step_time = step_time
        self.updated = False

    @property
    def min_dwell_time(self):
        """
        min_dwell_time : (s) min time between two successive steps
            if dwell time between two detected steps is shorter the second step
            will be ignored
        """
        return self._min_dwell_time

    @min_dwell_time.setter
    def min_dwell_time(self, min_dwell_time):
        self._min_dwell_time = min_dwell_time
        self.updated = False

    @property
    def steps(self):
        concatenated = np.r_[self._steps['pos'], self._steps['neg']]
        con_sort = np.sort(concatenated)
        return con_sort

    @property
    def pos_steps(self):
        return self._steps['pos']

    @property
    def neg_steps(self):
        return self._steps['neg']

    @property
    def plateau_segments(self):
        start = None
        stop = None
        if self.include_first_last:
            start = 0
            stop = self.max_stop
        return sn.idx_to_idx_segments(self.steps, start=start, stop=stop)

    @property
    def pos_plateau_segments(self):
        plateau_segments = self.plateau_segments
        idx = np.in1d(plateau_segments[:, 0], self.pos_steps)
        return plateau_segments[idx]

    @property
    def neg_plateau_segments(self):
        plateau_segments = self.plateau_segments
        idx = np.in1d(plateau_segments[:, 0], self.neg_steps)
        return plateau_segments[idx]


class CycleSectioner(SignalFeature):
    """
    A `CycleSectioner` analyses up to two signals for the excited signal ('x',
    'y'), the extrema ('minima', 'maxima'), the direction, i.e. positive
    or negative signal ('right', 'left'), and different cycles, i.e. rising
    ('rightstress', 'leftrelease') and falling ('rightrelease', 'leftstress')
    respectively.
    """
    def __init__(self, signals, resolution, compare_time=0.01,
                 threshold=0.005, min_duration=0.2, highest_frequency=1,
                 reduce_false_extrema=False, **kwargs):
        """
        Parameters
        ----------
        signals : 2D numpy.ndarray
            One or two signals to be analysed.
        resolution : float
            Number of datapoints/s of the resolution of the signals.
        compare_time : float, optional
            Time in s for size of moving window used to compare amplitude of
            signals to each other.
        threshold : float, optional
            Value, the amplitude of a signal is compared to, if only one signal
            is provided.
        min_duration : float, optional
            Minimum duration in s a preliminary detected excited section has to
            have to be finally accepted as one.
        highest_frequency : float, optional
            Highest expected frequency in Hz of occurence of extrema in the
            signals. If you set a too low frequency, extrema will not be
            detected. A frequency to high will result in to many falsly
            positive detected extrema. See also parameter
            `reduce_false_extrema`.
        reduce_false_extrema : bool
            Set to True, to reduce falsly positive detected extrema. As a side
            effect, the extrema within half a wavelength of the given highest
            frequency at the ends of the signal could be shifted or even not be
            detected. Make sure to properly set `highest_frequency`, first.
            Another option to reduce False positives is to smoothen the
            signals, before.

        Attributes
        ----------
        excited : dict
            Use 'x' and 'y' keys to get excited segments for the corresponding
            signal axis.
        sections : dict
            Use 'x' and 'y keys for corresponding axis. This will give you a
            list of dicts, containing the extrema indices with the keys
            'minima' and 'maxima' and the section segments with the keys
            'rightstress', 'rightrelease', 'leftstress', 'leftrelease',
            'right', 'left', 'rise', and 'fall'.
        compare_time
        threshold
        min_duration
        highest_frequency
        reduce_false_extrema
        """
        self._compare_time = compare_time
        self._threshold = threshold
        self._min_duration = min_duration
        self._highest_frequency = highest_frequency
        self._reduce_false_extrema = reduce_false_extrema

        super().__init__(signals=signals, resolution=resolution)

    def _update(self, signals, resolution):
        # excited sections for axes 'x' and 'y'
        self.excited = {}
        # sections (minima, maxima, left, right, stress, release, rise, fall)
        # of excited sections for axes 'x' and 'y'
        self.sections = {}
        self._set_sections(signals, resolution)

    def _set_sections(self, signals, resolution):
        # Set excited sections
        self.excited['x'], self.excited['y'] \
            = sn.get_excited_sections(signals, resolution,
                                      threshold=self.threshold,
                                      compare_time=self.compare_time,
                                      min_duration=self.min_duration)[0]

        # Set sections with information of left, right, release, and stress,
        s = {}
        if signals.ndim >= 2 and signals.shape[1] >= 2:
            s['x'] = signals[:, 0]
            s['y'] = signals[:, 1]
        else:
            if signals.ndim == 1:
                s['x'] = signals
            else:
                s['x'] = signals[:, 0]
            s['y'] = np.empty((0, 2), dtype=int)
        for axis in ['x', 'y']:
            self.sections[axis] = []
            for section in self.excited[axis]:
                try:
                    signal = s[axis][section[0]:section[1]]
                    minima, maxima = sn.get_extrema(signal, resolution,
                                                    self.highest_frequency,
                                                    self.reduce_false_extrema)
                    posrise, posfall, negfall, negrise, pos, neg, rise, fall \
                        = sn.get_sections(signal, minima, maxima)[0]
                    shift = section[0]
                    seg_desc = {'minima': minima + shift,
                                'maxima': maxima + shift,
                                'rightstress': posrise + shift,
                                'rightrelease': posfall + shift,
                                'leftstress': negfall + shift,
                                'leftrelease': negrise + shift,
                                'right': pos + shift,
                                'left': neg + shift,
                                'rise': rise + shift,
                                'fall': fall + shift}
                    self.sections[axis].append(seg_desc)
                except:
                    warnings.warn("Could not determine the minima and maxima "
                                  "of the excited position.")

    @property
    def compare_time(self):
        return self._compare_time

    @compare_time.setter
    def compare_time(self, compare_time):
        self._compare_time = compare_time
        self.updated = False

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, threshold):
        self._threshold = threshold
        self.updated = False

    @property
    def min_duration(self):
        return self._min_duration

    @min_duration.setter
    def min_duration(self, min_duration):
        self._min_duratione = min_duration
        self.updated = False

    @property
    def highest_frequency(self):
        return self._highest_frequency

    @highest_frequency.setter
    def highest_frequency(self, highest_frequency):
        self._highest_frequency = highest_frequency
        self.updated = False

    @property
    def reduce_false_extrema(self):
        return self._reduce_false_extrema

    @reduce_false_extrema.setter
    def reduce_false_extrema(self, reduce_false_extrema):
        self._reduce_false_extrema = reduce_false_extrema
        self.updated = False
