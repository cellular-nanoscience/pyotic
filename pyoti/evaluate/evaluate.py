# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:55:34 2016

@author: Tobias Jachowski
"""
import collections
import numpy as np
import matplotlib.pyplot as plt
import os

from . import signal as sn
from .. import gui
from .signalfeature import SignalFeature


class Evaluator(object):
    """
    An Evaluator extracts information from data provided by a Region and
    provides functions to further process and evaluate the data.

    The Evaluator holds references to one Region and one SignalFeature. The
    Region provides the data (amongst other specific details) to be evaluated.
    The SignalFeature analyses a signal from the data and identifies specific
    sections in the signal (like steps, plateaus, etc.). The Evaluator uses
    these identified sections to extract further information and processed data
    from the data provided by the Region.

    The Evaluator provides convenient functions to adjust the resolution of the
    (extracted) data and identified sections it returns.
    """
    def __init__(self, region=None, calibration=None, resolution=None,
                 filter_time=-1.0, sf_class=None, traces_sf=None,
                 resolution_sf=None, filter_time_sf=None, **kwargs):
        """
        Parameters
        ----------
        region
        calibration : pyoti.calibration.calibration.Calibration, optional
        resolution : float, optional
            Desired resolution of data (and sections) returned by this
            Evaluator in datapoints per s. The resolution is automatically
            limited to be between 1 and `region.samplingrate` datapoints per s.
            Additionally, resolution is limited to values of
            `self.samplingrate` / decimate, where decimate is a positive
            integer. The actual resolution can be checked with the attribute
            `self.resolution'.The resolution defaults to `region.samplingrate`.
        filter_time : float, optional
            Desired time containing datapoints to be used with a running
            filter. If `filter_time` == 0.0, averaging is switched off. If
            `filter_time` > 0.0, signal is filtered according to the given
            `filter_time`. If `filter_time` < 0.0 (default), signal is filtered
            according to `resolution` (`self.decimate`). The actual datapoints
            beeing used to filter the data can be checked with the property
            `self.window`.
        sf_class : pyoti.evaluate.signalfeature.SignalFeature
            The class of the SignalFeature to be used for signal analysis.
        traces_sf : str, int, list of str/int, slice
            The traces of the signal, the SignalFeature needs to analyze.
        resolution_sf : float, optional
            The resolution the SignalFeature should work with.
            Defaults to `resolution` (or `region.samplingrate`). Actual
            reslution can be checked with the attribute `self.resolution_sf`.
        filter_time_sf : float, optional
            The time for the moving filter, the signal for the SignalFeature
            should be filtered with. See parameter `filter_time` for details.
            Defaults to `filter_time`. The actual datapoints beeing used to
            filter the signal for the SignalFeature can be checked with the
            property `self.window_sf`.

        Attributes
        ----------
        traces_sf : str, int, list of str/int, slice
        region : pyoti.region.region.Region
        resolution : float
            Actual resolution of data (and sections) returned by this
            Evaluator in datapoints per second.
        moving_filter : str
            The moving filter to be applied for `self.get_data()`. Defaults to
            'mean'.
        filter_time : float
            Desired time containing datapoints to be used with a running
            filter. If `filter_time` == 0.0, averaging is switched off. If
            `filter_time` > 0.0, signal is filtered according to the given
            `filter_time`. If `filter_time` < 0.0, signal is filtered according
            to `resolution` (`self.decimate`). The actual datapoints beeing
            used to filter the data can be checked with the property
            `self.window`.
        window : int
            Actual datapoints beeing used to filter the data.
        resolution_sf : float
            Actual resolution the SignalFeature works with.
        moving_filter_sf : str
            The moving filter to be applied to the signal for the
            SignaleFeature. Defaults to 'mean'.
        filter_time_sf : float
            The time for the moving filter, the signal for the SignalFeature
            should be filtered with. See attribute `filter_time` for details.
            Defaults to `filter_time`. The actual datapoints beeing used to
            filter the signal for the SignalFeature can be checked with the
            property `self.window_sf`.
        window_sf : int
            Actual datapoints beeing used to filter the signal for the
            SignalFeature.
        timevector
        decimate
        decimate_sf
        calibration
        """
        # resolution for (almost) all properties, e.g. displacement, force, ...
        # self.update() once only, by setting self.region = region!
        # resolution 1.000 Hz, samplingrate 40.000 Hz -> 40 points

        self._resolution = resolution or region.samplingrate  # Hz
        self.moving_filter = 'mean'
        self._filter_time = filter_time  # s
        self._resolution_sf = resolution_sf or self._resolution
        self.moving_filter_sf = 'mean'
        self._filter_time_sf = filter_time_sf or filter_time

        self._sf_class = sf_class or SignalFeature
        self._sf = None

        self.traces_sf = traces_sf
        self.region = region  # implicitly calls self.update()

        self.calibration = calibration

        self.rfigure = None

    def get_data(self, traces=None, samples=None, copy=False, resolution=0,
                 moving_filter=None, filter_time=-1.0):
        """
        Get data of `self.region` with a defined resolution and smoothed by
        a moving filter.

        Parameters
        ----------
        traces : str or list of str, optional
            see parameters of method `pyoti.region.Region.get_data()`
        samples : int or list of int or slice, optional
            Index of samples that should be returned. The index is relative to
            the data returned by this method without any parameter given. If
            you want to have a defined resolution, via the paramters
            `resolution` or the attribute `self.resolution`, do only use int or
            slice and do set the step value of the slice to None or 1.
        copy: bool, optional
            see parameters of method `pyoti.region.Region.get_data()`
        resolution : int, optional
            Set the resolution of the returned data. If `resolution` <= 0
            (default), use `self.resolution` to set the resolution of the data.
        moving_filter : str, optional
            The moving filter to apply to the data. Can be one of 'mean' or
            'median'. Only applied, if `filter_time` != 0.0. Defaults to
            `self.moving_filter`.
        filter_time : float, optional
            Length of the window in s used to smooth the data with a moving
            filter. If `filter_time` is < 0.0 (default), use
            `self.filter_time` as the window length. If `filter_time` == 0.0,
            do not smooth the data.

        Returns
        -------
        2D numpy.ndarray
        """
        if resolution <= 0:
            decimate = self.decimate
        else:
            decimate = self._get_decimate(resolution)

        moving_filter = moving_filter or self.moving_filter

        if filter_time < 0.0:
            window = self.window
        elif filter_time > 0.0:
            window = max(1, int(filter_time * self.region.samplingrate))
        else:
            # filter_time == 0.0 # no averaging
            window = 1

        # Convert samples into the scope of self.region.get_data()
        samples = self.decimate_and_limit(samples)

        return self.region.get_data(traces=traces, samples=samples, copy=copy,
                                    moving_filter=moving_filter, window=window,
                                    decimate=decimate)

    def update(self):
        window = self.window_sf
        decimate = self.decimate_sf
        resolution = self.resolution_sf
        moving_filter = self.moving_filter_sf

        signals = self.region.get_data(self.traces_sf, copy=False,
                                       moving_filter=moving_filter,
                                       window=window, decimate=decimate)

        if self._sf is None:
            self._sf = self._sf_class(signals, resolution=resolution)
        else:
            self._sf.update(signals, resolution)

    def sections(self, slices=True, decimate=None, range_concat=False,
                 info=False, **kwargs):
        """
        Parameters
        ----------
        slices : bool, optional
            Return slices in a list (True) or a 2D numpy.ndarray with
            start/stop pairs in rows (False). Only evaluated, if `range_concat`
            is False.
        decimate : int, optional
            The step value for the slices returned. Only evaluated, if `slices`
            is True and `range_concat` is False.
        range_concat : bool, optional
            Return one 1D numpy.ndarray with indices, instead of individual
            slices or segments.
        info : bool, optional
            Additionally to the sections, return a 2D numpy.ndarray, containing
            infos about the individual sections (like which axis), if `info` is
            implemented in `self._sections()`. For details about the parmater
            see `self._sections()`.
        **kwargs, optional
            See method `self._sections()` for possible parameters.
        """
        # call function usually overwritten by subclasses
        sections = self._sections(info=info, **kwargs)
        if info and isinstance(sections, tuple):
            infos = sections[1]
            sections = sections[0]

        # Sort sections according to first index
        if sections.ndim == 1:
            sort = sections.argsort(axis=0)
        else:
            sort = sections.argsort(axis=0)[:, 0]
        sections = sections[sort]

        if info:
            infos = infos[sort]

        # undecimate and limit the sections
        sections = self.undecimate_and_limit(sections)

        # Convert segments (only, no idx!) to slices
        if slices and not range_concat and sections.ndim >= 2:
            sections = sn.idx_segments_to_slices(sections, decimate=decimate)

        # Convert segments (only, no idx!) to concatenated ranges
        if range_concat and sections.ndim >= 2:
            sections = sn.idx_segments_to_range(sections)

        if info:
            return sections, infos
        return sections

    def undecimate_and_limit(self, sections):
        # undecimate the sections
        undecimate = self.decimate_sf / self.decimate
        max_stop = int(np.ceil(self.region.datapoints / self.decimate))
        sections = np.round(sections * undecimate).astype(int)
        sections = sn.limit_segments(sections, min_start=0, max_stop=max_stop)
        return sections

    def decimate_and_limit(self, samples):
        if samples is None:
            return None

        decimate = self.decimate
        max_stop = self.region.datapoints

        # Convert different samples types into nunmpy.ndarray, which can be
        # checked by the function `sn.limit_segments`.
        if isinstance(samples, int):
            _samples = np.array([samples])
        elif isinstance(samples, slice):
            _samples = np.array([[samples.start or 0,
                                 samples.stop or self.datapoints]])
            step = samples.step
            if step is not None:
                step = np.round(samples.step * decimate).astype(int)
        else:
            _samples = samples
        _samples = np.round(_samples * decimate).astype(int)
        _samples = sn.limit_segments(_samples, min_start=0, max_stop=max_stop)

        if isinstance(samples, slice):
            samples = slice(_samples[0, 0], _samples[0, 1], step)
        else:
            samples = _samples

        return samples

    def _sections(self, **kwargs):
        sections = np.empty((0, 2), dtype=int)
        # create sections (segments or idx)
        return sections

    def init_rfig(self, legend=True, show=True, plot_params=None):
        if self.rfigure is not None:
            plt.close(self.rfigure)

        # Initialize proper plot parameters
        gui.set_plot_params(plot_params=plot_params)

        # create static figure
        self.rfigure = self._rfigure(legend=legend)

        if show:
            self.rfigure.show()

    def _rfigure(self, legend=True):
        figure = None
        # create the figure ...
        return figure

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, region):
        if region is None:
            raise TypeError("Missing required argument: 'region'")
        self._region = region
        self.update()

    @property
    def resolution(self):
        resolution = self.region.samplingrate / self.decimate
        return resolution

    @resolution.setter
    def resolution(self, resolution):
        self._resolution = resolution
        self.update()

    @property
    def filter_time(self):
        return self._filter_time

    @filter_time.setter
    def filter_time(self, filter_time):
        self._filter_time = filter_time
        self.update()

    @property
    def window(self):
        if self.filter_time < 0.0:
            window = self.decimate
        elif self.filter_time > 0.0:
            window = max(1, int(self.filter_time * self.region.samplingrate))
        else:  # filter_time == 0.0 # no averaging
            window = 1
        return window

    @property
    def resolution_sf(self):
        resolution_sf = self.region.samplingrate / self.decimate_sf
        return resolution_sf

    @resolution_sf.setter
    def resolution_sf(self, resolution_sf):
        self._resolution_sf = resolution_sf
        self.update()

    @property
    def filter_time_sf(self):
        return self._filter_time_sf

    @filter_time_sf.setter
    def filter_time_sf(self, filter_time_sf):
        self._filter_time_sf = filter_time_sf
        self.update()

    @property
    def window_sf(self):
        if self.filter_time_sf < 0.0:
            window_sf = self._get_decimate(self._resolution_sf)
        elif self.filter_time_sf > 0.0:
            window_sf = max(1, int(self.filter_time_sf
                                   * self.region.samplingrate))
        else:  # window_sf == 0.0 # no averaging
            window_sf = 1
        return window_sf

    @property
    def timevector(self):
        return self.region.timevector[::self.decimate]

    @property
    def decimate(self):
        return self._get_decimate(self._resolution)

    @property
    def datapoints(self):
        return int(np.ceil(self.region.datapoints / self.decimate))

    @property
    def decimate_sf(self):
        return self._get_decimate(self._resolution_sf)

    def _get_decimate(self, resolution):
        return max(1, int(np.round(self.region.samplingrate
                                   / max(1, resolution))))

    @property
    def calibration(self):
        if hasattr(self, '_calibration') and self._calibration is not None:
            return self._calibration
        else:
            return self.region.calibration

    @calibration.setter
    def calibration(self, calibration):
        self._calibration = calibration


def save_figures(figures, directory=None, file_prefix=None, file_suffix=None,
                 file_extension='.png', index_digits=3):
    """
    Save matplotlib figures in a given directory.

    The filenames of the figures will be a concatenation of the `file_prefix`,
    an index with `index_digits` digits, the `file_suffix` and the
    `file_extension`.

    Parameters
    ----------
    figures : Iterable of matplotlib figures or one figure
        A list, array, generator or other Iterable type of matplotlib figures.
        If figures is only one matplotlib figure, no index will be included in
        the filename of the figure.
    directory : str
        The directory, the figures should be saved in. The directory will be
        created, if it does not exist.
    file_prefix : str, optional
        A prefix every filename of the saved figures should include.
    file_suffix : str, optional
        A suffix every filename of the saved figures should include.
    files_extension : str, optional
        The file extension (and type) to be used to save the figures (default
        '.png').
    index_digits : int, optional
        Digits to be used for the index in the filename of the figures.
    """
    directory = directory or os.path.join(".", "results")
    file_prefix = file_prefix or ""
    file_suffix = file_suffix or ""

    # create results dir
    if not os.path.exists(directory):
        os.makedirs(directory)

    # If only one figure, save it without an index
    if not isinstance(figures, collections.Iterable):
        filename = "%s%s%s" % (file_prefix, file_suffix, file_extension)
        figures.savefig(os.path.join(directory, filename))

    for idx, fig in enumerate(figures):
        format_string = "".join(("%s%.", str(index_digits), "i%s%s"))
        filename = format_string % (file_prefix, idx, file_suffix,
                                    file_extension)
        fig.savefig(os.path.join(directory, filename))
