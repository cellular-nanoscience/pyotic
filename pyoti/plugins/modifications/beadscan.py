# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:31:32 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from pyoti.modification.modification import Modification, GraphicalMod
from pyoti import helpers as hp
from pyoti import traces as tc


class IBeadscan(GraphicalMod):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._ax = None

    def _figure(self):
        """
        Show the fitted model functions that are used to correct apply_data.
        """
        # create new figure and axes for adjusting the plateaus
        plot_params = {'lines.linewidth': 3.0}
        self._set_plot_params(plot_params)
        figure, ax_ = plt.subplots(1, 2, sharex=True, sharey=True)

        plottrace = self.modification.traces_apply
        ax = dict(list(zip(plottrace, ax_)))
        self._ax = ax

        # [ psdX, psdY ]
        traces = self.modification.traces_apply
        if 'psdZ' not in traces:
            traces.append('psdZ')
            traces.sort()

        # create lines to plot the corrections
        for trace in plottrace:
            ax[trace].plot([0], [0], c=tc.color(trace), marker='.', ls='')
            ax[trace].plot([0], [0], 'k-')
            ax[trace].plot([0], [0], 'y-', linewidth=3)
            ax[trace].legend(('fitdata', 'fit+extrapolation', 'fit'),
                             loc='best')
            ax[trace].set_xlabel(tc.label('psdZ'))
            ax[trace].set_ylabel(tc.label(trace))
            ax[trace].ticklabel_format(useOffset=False)
            ax[trace].grid(True)

        figure.suptitle("Crosstalk correction function for psdX and psdY from "
                        "psdZ.")

        return figure

    def _update_fig(self, **kwargs):
        plottrace = self.modification.traces_apply
        traces = self.modification.traces_apply
        traces.append('psdZ')
        Z = traces.index('psdZ')

        data_apply = self.modification._get_data_apply(
            traces=traces, decimate=True, copy=False)
        data_based = self.modification._get_data_based(
            traces=traces, decimate=True, copy=False)

        # Calculate the sorted indices for psdZ for plotting a "smooth" line
        # psdX, psdY, and psdZ sorted by increasing psdZ
        sortidx_apply = np.argsort(data_apply[:, Z])
        sortidx_based = np.argsort(data_based[:, Z])

        # get modification for traces psdX and psdY
        mod_index = hp.overlap_index(self.modification.traces_apply, plottrace)
        modification_apply \
            = self._modification(data_apply, self.samples_apply_decimated,
                                 traces, mod_index)
        modification_based \
            = self._modification(data_based, self.samples_based_decimated,
                                 traces, mod_index)

        ax = self._ax

        for trace in plottrace:
            lines = ax[trace].lines
            T = traces.index(trace)
            trace_idx_apply = self.modification.view_apply.traces_to_idx(trace)
            trace_idx_based = self.modification.view_based.traces_to_idx(trace)
            lines[0].set_data(data_based[:, Z], data_based[:, T])
            lines[1].set_data(data_apply[sortidx_apply, Z],
                              modification_apply[sortidx_apply,
                                                 trace_idx_apply])
            lines[2].set_data(data_based[sortidx_based, Z],
                              modification_based[sortidx_based,
                                                 trace_idx_based])

            # recompute ax.dataLim
            ax[trace].relim()
            # update ax.viewLim using new dataLim
            ax[trace].autoscale_view()

    def _modification(self, data, samples, data_traces, mod_index):
        """
        For compatibility reasons
        """
        psdZ_apply = self.modification._get_trace_data_apply(
            'psdZ', data, samples, data_traces, copy=False)

        # initialize the modification according to previously determined model
        # for X/Y
        mod_trace = hp.listify(mod_index)
        modification = np.zeros((len(data), len(mod_trace)))

        for i, trace in enumerate(mod_trace):
            modification[:, i] = self.modification._model[trace](psdZ_apply)

        return modification

    @property
    def samples_based_decimated(self):
        return slice(0, self.modification.view_based.datapoints,
                     self.modification.decimate)

    @property
    def samples_apply_decimated(self):
        return slice(0, self.modification.view_apply.datapoints,
                     self.modification.decimate)


class Beadscan(Modification):
    """
    Corrects for crosstalk from psdZ to psdX/psdY with a supplied beadscan
    in z axis.
    """
    GRAPHICALMOD = IBeadscan

    def __init__(self, **kwargs):
        traces_apply = ['psdX', 'psdY']
        super().__init__(datapoints=2500, traces_apply=traces_apply, **kwargs)

        # Define some attributes needed for the fitting of the beadscan
        # model for fitting psdX and psdY, usually UnivariatSpline
        self._model = {}

        # min and max value of psdZ used to determine the model
        self.minpsdZ = None
        self.maxpsdZ = None

    def _recalculate(self):
        # calculate data for fitting
        traces = self.traces_apply
        traces.append('psdZ')
        bin_means, bin_size = self.calculate_bin_means(traces=traces,
                                                       sorttrace=2)
        self.bin_means = bin_means

        # calculate the model for smoothing data
        self._calculate_model()

    def _print_info(self):
        for trace in self.traces_apply:
            print("    %s coeffs: %s"
                  % (tc.label(trace),
                     self._model[self.lia(trace)].get_coeffs()))

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        psdZ_apply = self._get_trace_data_apply('psdZ', data, samples,
                                                data_traces, copy=False)

        # apply modification and return modified data according to previously
        # determined model for X/Y
        data_trace = hp.listify(data_index)
        mod_trace = hp.listify(mod_index)

        for dtrace, mtrace in zip(data_trace, mod_trace):
            data[:, dtrace] -= self._model[mtrace](psdZ_apply)

        return data

    def _get_trace_data_apply(self, trace, data, samples, data_traces,
                              copy=True):
        # get psdZ and reduce dimensionality for later calculation of
        # modification
        if trace in data_traces:  # trace already in data
            index = data_traces.index(trace)
            trace_apply = data[:, index]
        else:  # psdZ needs to be fetched extra
            trace_apply = self._get_data_apply(samples=samples, traces=trace,
                                               copy=False)[:, 0]  # [ trace ]

        # check for extrapolation limit of the previously calculated model
        minpsdZ = trace_apply.min()
        maxpsdZ = trace_apply.max()
        if minpsdZ < self.minpsdZ:
            if self.minpsdZ - minpsdZ > 0.01:
                print("Beadscan: Need to extrapolate the data for crosstalk "
                      "psdZ -> psdX/Y min by %.5f V"
                      % (self.minpsdZ - minpsdZ))
        if maxpsdZ > self.maxpsdZ:
            if maxpsdZ - self.maxpsdZ > 0.01:
                print("Beadscan: need to extrapolate the data for crosstalk "
                      "psdZ -> psdX/Y max by %.5f V"
                      % (maxpsdZ - self.maxpsdZ))

        if copy:
            return trace_apply.copy()
        else:
            return trace_apply

    def _calculate_model(self):
        """
        Use beadscan x (psdX, psdY) and y (psdZ) data to create a
        one-dimensional smoothing model fit.
        Be aware of the indexing X, Y and Z (defined above) -> NOT psdX, psdY,
        psdZ
        """
        X = 0  # psdX
        Y = 1  # psdY
        Z = 2  # psdZ

        self.minpsdZ = self.bin_means[:, Z].min()
        self.maxpsdZ = self.bin_means[:, Z].max()
        self._model[X] = UnivariateSpline(self.bin_means[:, Z],
                                          self.bin_means[:, X])
        self._model[Y] = UnivariateSpline(self.bin_means[:, Z],
                                          self.bin_means[:, Y])
        self.set_changed()

    @property
    def bin_means(self):
        if not hasattr(self, '_v_bin_means'):
            return None
        else:
            return self._v_bin_means

    @bin_means.setter
    def bin_means(self, bin_means):
        self._v_bin_means = bin_means


# The following is only to update to database version 0.8.0
class GBeadscan(Beadscan):
    pass
