# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:29:19 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np

from . import signal as sn
from .evaluate import Evaluator
from .signalfeature import StepFinder


class Motion(Evaluator):
    """
    Calculates motion related parameters, like steps.

    To identify positive and negative steps, a Stepfinder is used.
    """
    def __init__(self, region=None, calibration=None, resolution=1000,
                 filter_time=-1.0, min_dwell_time=0.0, shift_time=0.0,
                 impose_dwell_time=0.0, resolution_sf=None,
                 filter_time_sf=None, **kwargs):
        """
        This __init__() constructor extends the superclass (Evaluator)
        constructor by setting the traces_sf to default to 'positionZ'.

        For extensive documentation of the Parameters and additional
        attributes, see `pyoti.evaluate.evaluate.Evaluator`.

        Parameters
        ----------
        region : pyoti.region.region.Region
        calibration : pyoti.calibration.calibration.Calibration
        resolution : float
        filter_time : float
        min_dwell_time : float
            Minimum time a step plateau needs to have.
        shift_time : float
            Time the start of the plateau is shifted.
        impose_dwell_time : float
            Adjust the stop index of the returned detected plateaus such a way,
            that they have a length of the given min_dwell_time.
        resolution_sf : float
        filter_time_sf : float
        **kwargs
            Extra arguments, e.g. traces_sf, in case the positionXY signals are
            called differently.

        Attributes
        ----------
        min_dwell_time
        min_plateau_length
        shift_time
        shift_length
        impose_dwell_time
        impose_plateau_length
        include_first_last_plateau
        pos_steps
        neg_steps
        steps
        pos_plateaus
        neg_plateaus
        plateaus
        """
        # resolution for (almost) all properties
        # resolution 1.000 Hz, samplingrate 40.000 Hz -> 40 points

        traces_sf = kwargs.pop('traces_sf', 'positionZ')

        super().__init__(region=region, calibration=calibration,
                         resolution=resolution, filter_time=filter_time,
                         sf_class=StepFinder, traces_sf=traces_sf,
                         resolution_sf=resolution_sf,
                         filter_time_sf=filter_time_sf)

        self._min_dwell_time = min_dwell_time  # s
        self._shift_time = shift_time  # s
        self._impose_dwell_time = impose_dwell_time  # s

    def _sections(self, plateaus=True, direction=None, **kwargs):
        if plateaus:
            if direction is 'pos':
                sections = self._sf.pos_plateau_segments
            elif direction is 'neg':
                sections = self._sf.neg_plateau_segments
            else:
                sections = self._sf.plateau_segments
            sections = self._shift_and_size_segments(sections)
        else:  # steps instead of plateaus
            if direction is 'pos':
                steps = self._sf.pos_steps
            elif direction is 'neg':
                steps = self._sf.neg_steps
            else:
                steps = self._sf.steps
            steps = steps + self.shift_length
            steps = sn.limit_segments(steps, min_start=0,
                                      max_stop=self.region.datapoints)
            sections = steps

        return sections

    def _shift_and_size_segments(self, segs):

        if self.impose_dwell_time > 0.0:
            # Select the plateaus having at least a size of min_dwell_time or
            # shift_time and impose_dwell_time
            size = max(self.min_plateau_length, self.shift_length
                       + self.impose_plateau_length)
            segs = sn.get_min_size_segments(segs, size)
            # Shift and size the plateaus and eventually crop the first/last
            # plateaus to min_start and max_stop
            segs = sn.shift_and_size_segments(segs, shift=self.shift_length,
                                              size=self.impose_plateau_length,
                                              min_start=0,
                                              max_stop=self.region.datapoints)
        else:
            # Select the plateaus having at least a size of min_dwell_time
            segs = sn.get_min_size_segments(segs, self.min_plateau_length)
            # Shift the plateaus and eventually crop the first/last plateaus to
            # min_start and max_stop
            segs = sn.shift_segments(segs, shift=self.shift_length,
                                     min_start=0,
                                     max_stop=self.region.datapoints)

        # After eventually having cropped the segments, reselect the plateaus
        # having at least a length of min_dwell_time
        segs = sn.get_min_size_segments(segs, self.min_plateau_length)

        return segs

    def _rfigure(self, legend=True):
        # Create figure and close old one
        figure = plt.figure()
        ax = figure.gca()

        ax.grid(True)
        ax.hold(True)

        # Get t, positionZ, psdXYZ and detected plateaus
        t = self.timevector
        positionZ = self.get_data(self.traces_sf)
        plateaus = self.plateaus

        # Plot the original signal
        ax.plot(t, positionZ, color='k', alpha=0.1)

        # Plot the detected plateaus (positionZ) and store psd data chunks and
        # heights
        for plat in plateaus:
            ax.plot(t[plat], positionZ[plat])

        ax.set_ylabel("Signal %s (um)" % (self.traces_sf))
        ax.set_xlabel('Time (s)')
        figure.suptitle('Automatically detected plateaus.')

        return figure

    @property
    def min_dwell_time(self):
        return self._min_dwell_time

    @min_dwell_time.setter
    def min_dwell_time(self, min_dwell_time):
        self._min_dwell_time = min_dwell_time
        self.update()

    @property
    def min_plateau_length(self):
        return int(np.round(self.min_dwell_time * self.resolution))

    @property
    def shift_time(self):
        return self._shift_time

    @shift_time.setter
    def shift_time(self, shift_time):
        self._shift_time = shift_time
        self.update()

    @property
    def shift_length(self):
        return int(np.round(self.shift_time * self.resolution))

    @property
    def impose_dwell_time(self):
        return self._impose_dwell_time

    @impose_dwell_time.setter
    def impose_dwell_time(self, impose_dwell_time):
        self._impose_dwell_time = impose_dwell_time
        self.update()

    @property
    def impose_plateau_length(self):
        return int(np.round(self.impose_dwell_time * self.resolution))

    @property
    def include_first_last_plateau(self):
        return self._sf.include_first_last

    @include_first_last_plateau.setter
    def include_first_last_plateau(self, include):
        self._sf.include_first_last = True
        self.update()

    @property
    def pos_steps(self):
        return self.sections(plateaus=False, direction='pos')

    @property
    def neg_steps(self):
        return self.sections(plateaus=False, direction='neg')

    @property
    def steps(self):
        return self.sections(plateaus=False)

    @property
    def pos_plateaus(self):
        return self.sections(direction='pos')

    @property
    def neg_plateaus(self):
        return self.sections(direction='neg')

    @property
    def plateaus(self):
        return self.sections()
