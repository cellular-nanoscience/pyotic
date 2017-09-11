# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 15:23:27 2017

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from pyoti.modification.modification import Modification, GraphicalMod
from pyoti import helpers as hp
from pyoti import traces as tc
from pyoti.evaluate import signal as sn
from pyoti.evaluate.tether import Tether


class IBaseline(GraphicalMod):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _figure(self):
        # Create figure and axes to plot
        plot_params = {'wspace': 0.09, 'hspace': 0.1}
        self._set_plot_params(plot_params)
        fig, ax = plt.subplots(2, 2)
        ax = ax.flatten()

        # create lines to plot the data
        for i, trace in enumerate(self.modification.traces_apply):
            # original data
            ax[i].plot([0], [0], 'c')
            # corrected data
            ax[i].plot([0], [0], 'm')
            # plot the model
            ax[i].plot([0], [0], 'r')
            # datapoints, which were used to fit the model
            ax[i].plot([0], [0], 'r.', ms=3.5)

            ax[i].legend(('original', 'corrected', 'fit', 'fit points'),
                         loc='best')
            ax[i].set_xlabel("Time in (s)")
            ax[i].set_ylabel(tc.label(trace))
            ax[i].grid(True)

        fig.suptitle("Original and corrected data and datapoints used.")

        return fig

    def _update_fig(self, **kwargs):
        # get the axes to plot
        ax = self.figure.axes

        # plot with 200 Hz
        decimate = int(np.round(self.modification.view_based.samplingrate
                                / 200))
        # average with 1000 Hz
        window = int(np.round(self.modification.view_based.samplingrate
                              / 1000))

        # Get the data to be plotted
        t = self.modification.view_based.timevector[::decimate]
        data = self.modification._get_data_based(
            traces=self.modification.traces_apply, copy=False,
            decimate=decimate, window=window)
        t_f, data_f = self.modification._calculate_fit_data()

        # Plot the data
        for i, trace in enumerate(self.modification.traces_apply):
            # original data
            ax[i].lines[0].set_data(t, data[:, i])
            # corrected data
            ax[i].lines[1].set_data(t, data[:, i]
                                    - self.modification._model[trace](t))
            # plot the model
            ax[i].lines[2].set_data(t, self.modification._model[trace](t))
            # datapoints, which were used to fit the model
            ax[i].lines[3].set_data(t_f, data_f[:, i])

            # recompute ax.dataLim
            ax[i].relim()
            # update ax.viewLim using new dataLim
            ax[i].autoscale_view()

        # plot the results of the tether object
        ax[3].clear()
        tether = self.modification._create_tether()
        tether._rfigure(ax=ax[3])


class Baseline(Modification):
    """
    Corrects the baseline of psdX, psdY and, psdZ.

    This Modification's `view_based` and `view_apply` need to have the same
    parent Region!
    """
    GRAPHICALMOD = IBaseline

    def __init__(self, db_update=False, **kwargs):
        traces_apply = ['psdX', 'psdY', 'psdZ']
        super().__init__(traces_apply=traces_apply, automatic_switch=True,
                         **kwargs)

        description = "Automatically detect baseline index"
        self.add_iattribute('auto_detect_idx', description=description,
                            value=True, unset_automatic=False)

        # binning time to calculate the models in s
        description = "Time window (s) used to calculate the smoothing model"
        self.add_iattribute('bin_time', description=description,
                            value=0.01, unset_automatic=False)

        description = "Traces to be corrected"
        self.add_iattribute('modify_traces', description=description,
                            value=self.traces_apply, unset_automatic=False)

        # Initialize the baseline index to meaningfull values. The indices
        # mark the datpoints that are used to calculate the models.
        self._baseline_idx = None
        if not db_update:
            self._calculate_baseline_idx()

        # smoothing models for self.traces_apply (usually UnivariatSpline),
        # which are used to calculate "fitted" values for the baseline,
        # which are subtracted as correction from the traces.
        self._model = {}

    def _recalculate(self):
        # Update the baseline_idx which is used to calculate the means
        if self.iattributes.auto_detect_idx:
            self._calculate_baseline_idx()

        # Calculate the bin means for the data to be used
        t, means = self._calculate_fit_data()

        # Calculate the model for modifying the data
        self._calculate_model(t, means)

    def _calculate_baseline_idx(self):
        tether = self._create_tether()
        _idx = tether.baseline_idx(**self.baseline_kwargs)
        self._baseline_idx = tether.decimate_and_limit(_idx)

    def _create_tether(self):
        # Get tether arguments
        tether_kwargs = self.tether_kwargs
        return Tether(**tether_kwargs)

    def _calculate_fit_data(self):
        """
        Calculate the means of the baseline and return the times, where the
        means were calculated, and the means itself.
        """
        # Get the data to calculate the means from
        data = self._get_data_based(traces=self.traces_apply, copy=False)

        # Get the indices of where the bin means should be calculated from,
        # and take the datapoints around these indices with a length according
        # to self.bin_time
        resolution = self.view_based.samplingrate
        bin_time = self.iattributes.bin_time
        window = int(np.round(resolution * bin_time / 2))
        starts = self.baseline_idx - window
        stops = self.baseline_idx + window + 1
        segs = np.c_[starts, stops]
        max_stop = self.view_based.datapoints
        segs = sn.limit_segments(segs, max_stop=max_stop)
        slices = sn.idx_segments_to_slices(segs)

        # Get the indices of the time points which are used as x values for the
        # model fitting
        t = self.view_based.timevector[self.baseline_idx]

        # Finale calculate the means
        means = hp.calculate_means(data, slices)

        return t, means

    def _calculate_model(self, t, means):
        """
        Fit a model (UnivariatSpline) to the given data.

        t : 1D numpy.ndarray
        means : 2D numpy.ndarray
            The y data to calculate the model. The shape of `means` has to be
            (`t.size`, `len(self.traces_apply)`)
        """
        for i, trace in enumerate(self.traces_apply):
            y = means[:, i]
            us = UnivariateSpline(t, y, s=0, ext=1)
            self._model[trace] = us

        self.set_changed()

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        # Get the shift of the time of the views the models were
        # calculated from (view_based) and the data is modified with
        # (view_apply). It is essential, that view_based and view_apply have
        # the same parent.
        shift = self.view_based.tmin - self.view_apply.tmin

        # Get the timevector for the requested samples and correct it with the
        # timeshift. It is used by the UnivariatSpline
        t = self.view_apply.timevector[samples]
        t -= shift

        # Get the indices of the data traces to be modified and the
        # names of the traces to be modified
        data_indices = hp.listify(data_index)
        mod_traces = np.array(self.traces_apply)[mod_index]

        # Modify the data with the UnivariatSpline
        for didx, mtrace in zip(data_indices, mod_traces):
            # Check if trace should be modified
            if mtrace in self.iattributes.modify_traces:
                data[:, didx] -= self._model[mtrace](t)

        return data

    @property
    def baseline_idx(self):
        return self._baseline_idx

    @baseline_idx.setter
    def baseline_idx(self, baseline_idx):
        self._baseline_idx = np.unique(baseline_idx)
        self.iattributes.auto_detect_idx = False
        self.set_changed()

    @property
    def tether_kwargs(self):
        if not hasattr(self, '_tether_kwargs') or self._tether_kwargs is None:
            tether_kwargs = {}
        else:
            tether_kwargs = self._tether_kwargs.copy()
        # Default argument `region` to view_based
        tether_kwargs['region'] = tether_kwargs.get('region', self.view_based)
        return tether_kwargs

    @tether_kwargs.setter
    def tether_kwargs(self, tether_kwargs):
        self._tether_kwargs = tether_kwargs
        self.set_changed()

    @property
    def baseline_kwargs(self):
        if not hasattr(self, '_baseline_kwargs') \
                or self._baseline_kwargs is None:
            return {}
        baseline_kwargs = self._baseline_kwargs.copy()
        return baseline_kwargs

    @baseline_kwargs.setter
    def baseline_kwargs(self, baseline_kwargs):
        self._baseline_kwargs = baseline_kwargs
        self.set_changed()


# The following is only to update to database version 0.8.0
class GBaseline(Baseline):
    pass
