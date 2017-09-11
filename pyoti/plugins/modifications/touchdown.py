# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:28:51 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor
from scipy.optimize import fsolve

from pyoti.modification.modification import Modification, GraphicalMod
from pyoti import traces as tc


class ITouchdown(GraphicalMod):
    """
    Subclass of Touchdown that provides graphical interfaces to adjust the fit
    parameters.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Define attributes needed for the graphical fitting of the touchdown
        # a cursor on the axes to show upper fitting border and left/right
        # (middle) border
        self._cursor = None
        self._lines = {}
        self._ax = None
        self.supertitle = None

    def _figure(self):
        """
        Initialize and show an interactive plot for fitting the touchdown.
        Adjust the fit parameters interactively by displaying the fitted
        touchdown. The plot of the fitting is stored in self.figure.
        """
        # create new figure and axes for fitting the touchdown
        figure, ax = plt.subplots(1, 1, sharex=True, sharey=True)
        self._ax = ax

        # create lines to plot the data for the fit and for displaying fitting
        # result
        self._lines['orig_data'] = ax.plot([0], [0], tc.color('psdZ'))[0]
        self._lines['binned_data'] = ax.plot([0], [0], 'b.')[0]

        # cursor to display the "boundaries" for fitdata
        self._cursor = Cursor(ax, useblit=True, color='red', lw=1)

        # initialize all lines for the plot that represent the fit paramters
        # extrapolation of fit left
        self._lines['ext_left'] = ax.plot([0, 0], [0, 0], 'c', lw=1)[0]
        # extrapolation of fit right
        self._lines['ext_right'] = ax.plot([0, 0], [0, 0], 'c', lw=1)[0]
        # fit left
        self._lines['fit_left'] = ax.plot([0, 0], [0, 0], 'r', lw=3)[0]
        # fit right
        self._lines['fit_right'] = ax.plot([0, 0], [0, 0], 'g', lw=3)[0]
        # touchdown position
        self._lines['touchdown'] = ax.axvline(x=0.0, color='y', lw=2)
        # touchdown of bead, where positionZ equal to 0.0
        self._lines['touchdown_calib'] = ax.axvline(x=0.0, color='k', lw=2,
                                                    dashes=(10, 10))

        ax.set_xlabel(tc.label('positionZ'))
        ax.set_ylabel(tc.label('psdZ'))
        ax.ticklabel_format(useOffset=False)
        ax.grid(True)

        self.supertitle = figure.suptitle("Fit touchdown determination (move "
                                          "to upper and right/left limit, "
                                          "left/right click, repeat ...)")

        # method to update plot and perform the fit
        figure.canvas.mpl_connect('button_release_event',
                                  self._update_touchdown)

        return figure

    def _update_fig(self, data_changed=True, **kwargs):
        """
        Updates the interactive fitting plot: Set values of fit parameters,
        touchdown position and fitted lines
        """
        if data_changed:
            # data, [:,0] is positionZ, [:,1] is psdZ
            data = self.modification._get_data_based(
                traces=['positionZ', 'psdZ'], decimate=True, copy=False)
            self.modification.validate_fit_params(data=data)
            bin_means = self.modification.bin_means
            minx = data[:, 0].min() - 0.05
            maxx = data[:, 0].max() + 0.05
            miny = data[:, 1].min() - 0.05
            maxy = data[:, 1].max() + 0.05
            self._ax.set_xlim([minx, maxx])
            self._ax.set_ylim([miny, maxy])

            self._lines['orig_data'].set_data(data[:, 0], data[:, 1])
            self._lines['binned_data'].set_data(bin_means[:, 0],
                                                bin_means[:, 1])

        # get left and right datapoint means according to fit borders for
        # fitting
        bin_means = self.modification.bin_means
        left_means = self.modification.left_means
        right_means = self.modification.right_means

        self._lines['ext_left'].set_data(bin_means[:, 0], self.get_pv_left(
                                         bin_means[:, 0]))
        self._lines['ext_right'].set_data(bin_means[:, 0], self.get_pv_right(
                                          bin_means[:, 0]))
        self._lines['fit_left'].set_data(left_means[:, 0], self.get_pv_left(
                                         left_means[:, 0]))
        self._lines['fit_right'].set_data(right_means[:, 0], self.get_pv_right(
                                          right_means[:, 0]))
        self._lines['touchdown'].set_xdata([self.modification.touchdown,
                                            self.modification.touchdown])
        self._lines['touchdown_calib'].set_xdata(
            [self.modification._calibration.touchdown,
             self.modification._calibration.touchdown])

    def _pre_close_fig(self):
        # Store touchdown fit plot for documentation
        self.supertitle.set_text("Touchdown determined by polynomial fits")
        self._cursor.visible = False
        # update visibility status of self._cursor before removing reference
        # to it
        self._figure_canvas_draw()

        # Release unused memory
        self._cursor = None
        self._ax = None
        self._lines.clear()
        self.supertitle = None

    def _update_touchdown(self, event):
        """
        Gets called by events of the canvas object of the figure
        """
        self._set_fit_params(event)
        touchdown = self.modification.fit_touchdown()
        self.modification._set_touchdown_leave_auto(touchdown)
        self.update_fig(data_changed=False)

    def _set_fit_params(self, event):
        """
        Set the parameters for the fitting routine.
        """
        if event.button == 1:  # left mouse button -> adjust left fitdata
            # left border for fitdata
            self.modification.left = self._cursor.linev.get_xdata()[0]
            # max value (psdZ) for fitdata
            self.modification.left_upper = self._cursor.lineh.get_ydata()[0]
        else:  # right (or middle) mouse button -> adjust right fitdata
            # right border for fitdata
            self.modification.right = self._cursor.linev.get_xdata()[0]
            # max value (psdZ) for fitdata
            self.modification.right_upper = self._cursor.lineh.get_ydata()[0]

    def get_pv_left(self, x):
        return np.polyval(self.modification.pf_left, x)

    def get_pv_right(self, x):
        return np.polyval(self.modification.pf_right, x)


class Touchdown(Modification):
    """
    Determine the approximate touchdown
    The displacement sensitivity and stiffness of a calibration are always
    relative to and given absolute at the touchdown of the bead on the surface
    The positionZ is touchdown corrected by this modification.
    """
    GRAPHICALMOD = ITouchdown

    def __init__(self, fit_touchdown=True, **kwargs):
        traces_apply = ['positionZ']
        # the touchdown position determined by calibration and fitting
        super().__init__(automatic_switch=True, datapoints=2500,
                         traces_apply=traces_apply, **kwargs)

        # register a widget, giving a key, a function called upon change
        self.add_iattribute('touchdown', description="Touchdown (V)",
                            value=0.0)

        # Initially set touchdown to touchdown of calibration, if touchdown
        # should not be automatically fitted. This will automatically switch of
        # the automatic fitting.
        if not fit_touchdown:
            self.touchdown = self._calibration.touchdown

        # register a button to set the touchdown from calibration.touchdown
        self.add_iattribute('set_touchdown_calib',
                            description="Set touchdown from calibration",
                            callback_functions=[
                                self._set_touchdown_from_calibration,
                                self.evaluate])

        # Define attributes needed for the fitting of the touchdown
        self.left_upper = None  # max value (psdZ) for fitdata
        self.right_upper = None  # max value (psdZ) for fitdata
        self.left = None  # right/left border for fitdata
        self.right = None  # right/left border for fitdata

    def _recalculate(self):
        # calculate data for fitting
        self.validate_fit_params()

        # fit the touchdown
        touchdown = self.fit_touchdown()
        self._set_touchdown_leave_auto(touchdown)

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        data[:, data_index] -= self.touchdown
        return data

    def validate_fit_params(self, data=None):

        # self.init_fit_parameters()
        if data is None:
            data = self._get_data_based(traces=['positionZ', 'psdZ'],
                                        decimate=True, copy=False)

        positionzmin = data[:, 0].min()
        positionzmax = data[:, 0].max()
        middle = positionzmin + (positionzmax - positionzmin) / 2
        upper = data[:, 1].max()

        # First initialisation of fit boundaries
        if self.left is None or self.right_upper is None:
            self.left_upper = upper   # max value (psdZ) for fitdata
            self.left = middle        # right/left border for fitdata
        if self.right is None or self.right_upper is None:
            self.right = middle       # right/left border for fitdata
            self.right_upper = upper  # max value (psdZ) for fitdata

        bin_means, bin_size = self.calculate_bin_means(data=data)
        self.bin_means = bin_means

        if not np.any(self.left_means_idx):
            # first, change only the upper boundary and leave the left one
            # usually solves problems, when psdZ was changed by an offset
            self.left_upper = upper
            if not np.any(self.left_means_idx):
                # still no means, change left boundary, too
                self.left = middle
        if not np.any(self.right_means_idx):
            # first, change only the upper boundary and leave the right one
            # usually solves problems, when psdZ was changed by an offset
            self.right_upper = upper
            if not np.any(self.right_means_idx):
                # still no means, change right boundary, too
                self.right = middle

    def fit_touchdown(self):
        """
        Performs the fit.
        In essence, it fits 2 polynomials of 2nd degree to the left (bead
        on surface) and to the right (free bead) portion of the psdZ data
        The intercept of these two polynomials determine the touchdown of the
        bead with the surface.
        """
        if self.left_means.shape[0] > 0 and self.right_means.shape[0] > 0:
            def f_left(x):
                return np.polyval(self.pf_left, x)

            def f_right(x):
                return np.polyval(self.pf_right, x)

            def findIntersection(fun0, fun1, x0):
                # one could also solve formula for the intersection
                # analytically, but fsolve is more flexible
                return fsolve(lambda x: fun0(x) - fun1(x), x0)
            middle = self.left + (self.right - self.left)/2
            touchdown = findIntersection(f_left, f_right, middle)[0]

            return touchdown

    def _set_touchdown_leave_auto(self, touchdown):
        self.iattributes.set_value('touchdown', touchdown,
                                   leave_automatic=True)

    def _set_touchdown_from_calibration(self):
        self.touchdown = self._calibration.touchdown

    @property
    def touchdown(self):
        return self.iattributes.touchdown

    @touchdown.setter
    def touchdown(self, touchdown):
        self.iattributes.touchdown = touchdown

    @property
    def _calibration(self):
        # calibration, which is used to:
        # - set touchdown from calibration
        # - visualize calibration.touchdown in GTouchdown
        return self.view_apply.calibration

    @property
    def pf_left(self):
        # get left datapoint means according to fit borders for fitting
        left_means = self.left_means
        # polynomial 2nd order order for negative positionZ
        return np.polyfit(left_means[:, 0], left_means[:, 1], 2)

    @property
    def pf_right(self):
        # get right datapoint means according to fit borders for fitting
        right_means = self.right_means
        # polynomial 2nd order for positive positionZ
        return np.polyfit(right_means[:, 0], right_means[:, 1], 2)

    @property
    def left_means(self):
        return self.bin_means[self.left_means_idx]

    @property
    def right_means(self):
        return self.bin_means[self.right_means_idx]

    @property
    def left_means_idx(self):
        x_check = self.bin_means[:, 0] < self.left
        y_check = self.bin_means[:, 1] < self.left_upper
        return np.logical_and(x_check, y_check)

    @property
    def right_means_idx(self):
        x_check = self.bin_means[:, 0] > self.right
        y_check = self.bin_means[:, 1] < self.right_upper
        return np.logical_and(x_check, y_check)

    @property
    def bin_means(self):
        if not hasattr(self, '_v_bin_means'):
            return None
        else:
            return self._v_bin_means

    @bin_means.setter
    def bin_means(self, bin_means):
        self._v_bin_means = bin_means
