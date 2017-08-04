# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:09:09 2016

@author: Tobias Jachowski
"""
import collections
import numpy as np
try:
    import pandas as pd
    __pd__ = True
except ImportError:
    __pd__ = False
from abc import ABCMeta, abstractmethod, abstractproperty

from .. import helpers as hp
from .. import traces as tc
from ..graph import GraphMember
from ..evaluate import signal as sn


class Region(GraphMember, metaclass=ABCMeta):
    """
    A Region provides common methods, to define and reference regions on given
    data, and methods to access the data in the defined regions.
    A Region of the data to refer to can be choosen by setter methods
    (start/stop or tmin/tmax). The region can be defined as a timespan
    (seconds) or an indexspan (indexes of a numpy array).
    """

    def __init__(self, caching=False, **kwargs):
        """
        Parameters
        ----------
        caching : bool
        **kwargs

        Attributes
        ----------
        start
        stop
        tmin
        tmax
        timespan
        samplingrate
        datapoints
        seconds
        indexspan
        indexvector
        timevector
        num_traces
        caching
        traces
        data
        """
        super().__init__(**kwargs)
        self.caching = caching

    @abstractproperty
    def start(self):
        pass

    @abstractproperty
    def stop(self):
        pass

    @start.setter
    def start(self, start):
        pass

    @stop.setter
    def stop(self, stop):
        pass

    @property
    def tmin(self):
        return self.start / self.samplingrate

    @property
    def tmax(self):
        return (self.stop - 1) / self.samplingrate

    @tmin.setter
    def tmin(self, tmin):
        self.start = int(np.around(tmin * self.samplingrate))

    @tmax.setter
    def tmax(self, tmax):
        self.stop = int(np.around(tmax * self.samplingrate)) + 1

    @property
    def timespan(self):
        return (self.tmin, self.tmax)

    @timespan.setter
    def timespan(self, timespan):
        self.tmin, self.tmax = timespan

    def set_timespan(self, tmin, tmax):
        self.tmin = tmin
        self.tmax = tmax

    @abstractproperty
    def samplingrate(self):
        # points per second
        pass

    @property
    def resolution(self):
        return self.samplingrate

    @property
    def datapoints(self):
        return self.stop - self.start

    @property
    def seconds(self):
        return self.datapoints / self.samplingrate

    @property
    def indexspan(self):
        # return slice(self.start, self.stop)
        return slice(0, self.datapoints, 1)

    @property
    def indexvector(self):
        # return np.arange(self.start, self.stop)
        return np.arange(self.datapoints)

    @property
    def timevector(self):
        return self.indexvector / self.samplingrate

    @abstractproperty
    def num_traces(self):
        return self.data.shape[0]

    @property
    def caching(self):
        return self._caching

    @caching.setter
    def caching(self, caching):
        """
        If memory consumption is no concern, use caching to speed up data
        retrieval.
        """
        self._caching = caching
        if not caching:
            # delete cache to free memory, ZODB volatile
            self._v_data_cached = None

    def get_data(self, traces=None, samples=None, moving_filter='mean',
                 window=1, decimate=1, copy=True, pandas=False):
        """
        Returns data of this region.

        Parameters
        ----------
        traces : str or list of str, optional
        samples : int, list of int, or slice, optional
            Samples start/stop has to be in between 0 and self.datapoints. If
            start, stop or any other index is negative it will be converted to
            a positive index, referring to self.datapoints.
        moving_filter : str, optional
            The moving filter to apply to the data. Can be one of 'mean' or
            'median'. Only applied, if `window` > 1.
        window : int, optional
            The window size of the moving filter.
        decimate : int, optional
            Decimate the data by given value. If `samples` is of type slice and
            `samples.step` is explicitly set (i.e. not None), the step value
            takes precedence over `decimate`.
        copy: bool, optional
            Only if `self.caching` is True, `copy` has the effect of either
            returning a copy of the cached data (`copy`=True) or returning a
            reference to the cached data array.
            If `self.caching` is False, return the uncached data. Then, `copy`
            has no effect at all.
        pandas : bool, optional
            Return a pandas.DataFrame with proper indexes and colunmns.

        Returns
        -------
        2D numpy.ndarray or pandas.DataFrame
            If self.caching is True, return cached data.
        """
        # Set the traces and the samples to proper default values and/or
        # convert (normalize) them to proper formats
        traces_idx = self.traces_to_idx(traces)
        samples = self.samples_idx(samples, decimate)

        if window <= 1:
            # Return unfiltered data
            data = self._get_data(samples, traces_idx, copy)
        else:
            # Moving filter should be applied, filter first, then decimate.
            # Get data, disregarding decimating factor, but considering traces
            if isinstance(samples, slice):
                samples_start = samples.start
                samples_stop = samples.stop
            else:  # samples is numpy.ndarray
                samples_start = samples.min()
                samples_stop = samples.max()
            # Get requested and preceding and succeeding data, to avoid
            # boundary issues of the filter
            filter_start = max(0, samples_start - window)
            filter_stop = min(samples_stop + window, self.datapoints)
            filter_samples = slice(filter_start, filter_stop, 1)
            data = self._get_data(filter_samples, traces_idx, copy=False)

            # Filter the data
            data = hp.moving_filter(data, window, moving_filter=moving_filter)

            # Return requested samples, considering decimating factor
            if isinstance(samples, slice):
                start = samples_start - filter_start
                stop = samples_stop - filter_start
                samples = slice(start, stop, samples.step)
            else:  # samples is numpy.ndarray
                samples = samples - filter_start
            data = data[samples]

        if pandas and __pd__:
            return pd.DataFrame(data, index=self.timevector[samples],
                                columns=self.idx_to_traces(traces_idx))

        return data

    def _get_data(self, samples, traces_idx, copy=True):
        """
        Returns the data for the given samples and traces_idx.
        """
        # Retrieve cached data
        if self.caching:
            # If the cache is outdated or not yet created, do so
            self.update_cache()

            # Check, if traces and/or samples are instances of numpy.ndarray
            traces_is_ndarray = isinstance(traces_idx, np.ndarray)
            samples_is_ndarray = isinstance(samples, np.ndarray)

            # If at least one of the indices is an ndarray, advanced indexing
            # is triggered.
            if traces_is_ndarray or samples_is_ndarray:
                # If both, traces and samples are an ndarray, indexing two
                # dimensions in a single step would result in an error in
                # numpy. Therefore, separate the indexing into two steps.
                if traces_is_ndarray and samples_is_ndarray:
                    # First, select the samples.
                    data_ = self._v_data_cached[samples]

                    # Second, select the traces.
                    return data_[:, traces_idx]

                # Advanced indexing automatically returns a copy of the data.
                # There is no need to extra copy via ndarray.copy(), which
                # would cost some extra execution time.
                copy = False

            if copy:
                # Per default return a copy of the data cache to protect it
                # from unwanted in place modifications
                return self._v_data_cached[samples, traces_idx].copy()
            else:
                # Only, if explicitly asked for "copy=False", return a
                # direct reference to the cached data
                return self._v_data_cached[samples, traces_idx]

        # No caching enabled, return uncached data
        return self._get_data_uncached(samples, traces_idx, copy)

    def update_cache(self, force=False):
        """
        Create and update cached data.

        Parameters
        ----------
        force : bool
            Recalculate, even if self._v_data_cached is up to date.
        """
        # Check for whether an update of the cache is needed or not.
        if self.caching and (not hasattr(self, '_v_data_cached')
                             or self._v_data_cached is None
                             or force):
            # An update is needed, calculate data for self.indexspan and
            # all traces_idx and store it in the cache. ZODB volatile.
            self._v_data_cached = self._get_data_uncached(self.indexspan,
                                                          self.traces_to_idx(),
                                                          copy=True)

    @abstractmethod
    def _get_data_uncached(self, samples, traces_idx, copy=True):
        pass

    def member_changed(self, ancestor=True, calledfromself=False,
                       leave_cache=False, **kwargs):
        # If self triggered a change, delete cache and trigger an update of the
        # cache. `leave_cache` can prevent deleting of the cache.
        if calledfromself and not leave_cache:
            self._v_data_cached = None  # ZODB volatile
        # If an ancestor triggered a change, delete cache and trigger an update
        # of the cache. A triggered change of descendants is ignored.
        if not calledfromself and ancestor:
            self._v_data_cached = None  # ZODB volatile

        # Call method of superclass `GraphMember`
        super().member_changed(ancestor=ancestor,
                               calledfromself=calledfromself, **kwargs)

    def traces_available(self, traces=None):
        traces_idx = self.traces_to_idx(traces)
        return self.idx_to_traces(traces_idx)

    def traces_to_idx(self, traces=None):
        """
        return index/slice of trace/s
        if trace is None, default to all traces
        """
        # Default to all traces
        if traces is None:
            traces_idx = slice(0, self.num_traces, 1)
        # If traces is a slice, correct for negative and None values
        elif isinstance(traces, slice):
            # TODO: if step is negative, indices() evaluates stop to
            # -1, but it should be -1 -self.num_traces!
            start, stop, step = traces.indices(self.num_traces)
            traces_idx = slice(start, stop, step)
        else:
            traces_idx = []
            for trace in tc.normalize(traces):
                if trace in self.traces:
                    traces_idx.append(self.traces.index(trace))
                elif isinstance(trace, int) and trace < self.num_traces:
                    traces_idx.append(trace)

        return hp.slicify(traces_idx, length=self.num_traces)

    def idx_to_traces(self, traces_idx=None):
        if isinstance(traces_idx, slice):
            traces_idx = hp.listify(traces_idx)
        if traces_idx is None:
            traces = self.traces
        elif isinstance(traces_idx, collections.Iterable):
            traces = [self.traces[trace_idx]
                      for trace_idx in traces_idx
                      if trace_idx < self.num_traces]
        return traces

    def samples_idx(self, samples=None, decimate=1):
        """
        Returns a slice or an np.array.

        Parameters
        ----------
        samples : int or list of int or slice, optional
            If samples is None, return the slice spanning the whole Region:
            slice(0, `self.datapoints`, `decimate`).
            If `samples` is an instance of int, slice, list or np.array,
            negative indices (i_neg) are converted into positive ones (i_pos):
            i_pos = i_neg + `self.datapoints`
        decimate : int, optional
            The decimate value of the returned slice, or the value the index
            array is decimated with.
            If `samples` is of type slice and `samples.step` is not explicitly
            set (i.e. None), `decimate` takes precedence over the step value.
        """
        if samples is None:
            samples = slice(0, self.datapoints, decimate or 1)
        elif isinstance(samples, int):
            if samples < 0:
                # correct for negative indices in original samples
                # referring to self.datapoints
                samples += self.datapoints
            samples = slice(samples, samples + 1, decimate or 1)
        elif isinstance(samples, slice):
            # Correct samples for negative indices and None
            start, stop, step = samples.indices(self.datapoints)

            # If samples.step is not explicitly set (i.e. None), decimate takes
            # precedence of the step information (1 per default)
            step = samples.step or decimate or 1

            samples = slice(start, stop, step)
        else:
            # samples either a list or an np.array
            # convert list to array and apply decimate (1 per default)
            if isinstance(samples, list):
                samples = np.array(samples)[::decimate or 1]
            # correct for negative indices in original samples
            # referring to self.datapoints
            samples[samples < 0] += self.datapoints

        return hp.slicify(samples)

    def _excited(self, traces=None, index=True):
        """
        For backward compatibility in pyotc.

        Parameters
        ----------
        traces : list of str
            Two traces to be compared. Defaults to ['positionX', 'positionY'].

        Returns
        -------
        trace : str
            The excited trace.
        """
        traces = traces or ['positionX', 'positionY']
        data = self.get_data(traces=traces, copy=False)
        ex = sn.get_excited_signal(data)
        if index:
            return ex
        return traces[ex]

    @abstractproperty
    def traces(self):
        return

    @property
    def data(self):
        # return all data
        return self.get_data()

    def __getattr__(self, name):
        """
        Allow attributes to be used as trace selections for get_data
        """
        if name in self.traces or name in tc:
            # name is directly known by Region or
            # name is probably an alias, a shorthand notation, or a combination
            return self.get_data(traces=name)
        else:
            raise AttributeError(name)
