# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:07:32 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np
import time
from ipywidgets import interact, IntSlider
from IPython.core.display import Image, display
from matplotlib.widgets import SpanSelector

from . import helpers as hp
from . import traces as tc


class GRS(object):
    """
    The GRS (Graphical Region Selector) creates a user interface to show data
    and select a range of the data.

    Attributes
    ----------
    ifigure : matplotlib.figure.Figure
        The current figure.
    displayrate : float
        Resolution of the plotted data in Hz.
    """

    def __init__(self, **kwargs):
        self._dpi = 100.0  # pixel per inch
        self._width = 800  # pixel per column
        self._height = 150  # pixel per row
        self._left = 65  # pixel for left
        self._right = 20  # pixel for right
        self._bottom = 55  # pixel for bottom
        self._top = 55  # pixel for top
        self._hspace = 18  # pixel between rows
        self._wspace = 18  # pixel between columns
        self.ifigure = None
        # _axspan is used by adjust_region AND _set_time_span, list of
        # indicators, which span is going to be used for every section
        self._axspan = {}
        # self._ax_OK = None
        self._tmin = 0
        self._tmax = 0
        self._minspan = 0.1  # s
        self._pre_title = "Select timespan (leftclick, draw, release)"
        self.displayrate = 1000.0  # display max 1000.0 datapoints per second
        self._onselect_callbacks = []
        self._spanselector = {}
        self._axspan = {}
        self.select = True

    def add_onselect_callback(self, function=None):
        """
        Register a function as a callback, when the selection of the data is
        changed.

        Parameters
        ----------
        function : function or method
            The function needs to accept two parameters tmin and tmax

        Returns
        -------
        bool
            True, if function was registered or is already registered as a
            callback function. Otherwise False.
        """
        if function is not None and function not in self._onselect_callbacks:
            self._onselect_callbacks.append(function)
            return True
        return False

    def set_onselect_callback(self, function=None):
        """
        Set a function as a callback, when the selection of the data is
        changed and deregister all other callback functions.

        Parameters
        ----------
        function : function or method
            The function needs to accept two parameters tmin and tmax.

        Returns
        -------
        bool
            True, if function was set as the callback function. Otherwise
            False.
        """
        if function is not None:
            self._onselect_callbacks.clear()
            self._onselect_callbacks.append(function)
            return True
        return False

    def _ifigure(self, timevector, data, samplingrate, traces, tmin=None,
                 tmax=None, xlim=None, description=None):
        """
        Create a figure based on the given timevector and data that lets the
        user choose a timespan (tmin, tmax):
          - Create figure with axes for all traces
          - Add SpanSelector for selecting timespan and set initial tmin and
            tmax according to given region
          - return figure

        Parameters
        ----------
        timevector : np.array
            1D array with the timevector to plot
        data : np.array
            The data array containing the traces to plot
        samplingrate : float
            The samplingrate of the data in Hz
        traces : list of str
            List with names of traces to use for data, COLOR and NAME
        tmin : float, optional
            tmin for the span selector
        tmax : float, optional
            tmax for the span selector
        xlim : tuple, optional
            tuple with (tmin, tmax) to set the range of the data to plot
        description : str, optional
            str describing the purpose of the region to select. Is used as a
            heading

        Returns
        -------
        matplotlib.figure.Figure
        """

        traces = tc.normalize(traces)

        rows = len(traces)
        columns = 1

        # Create figure and axes
        ifigure, axes_ = plt.subplots(rows, columns, sharex=True)
        if rows * columns == 1:
            axes_ = np.array([axes_])
        ax = dict(list(zip(traces, axes_)))

        if tmin is None:
            tmin = timevector[0]
        if tmax is None:
            tmax = timevector[-1]
        self._tmin = tmin
        self._tmax = tmax

        if xlim is None:
            xlim = (timevector[0], timevector[-1])

        decimate = max(1, int(np.round(samplingrate / self.displayrate)))

        # plot traces
        for i, trace in enumerate(traces):
            ax[trace].plot(timevector[::decimate], data[:, i][::decimate],
                           tc.color(trace))
            ax[trace].set_ylabel(tc.label(trace))
            ax[trace].ticklabel_format(useOffset=False)
            ax[trace].set_xlim(*xlim)
            ax[trace].grid(True)

            if self.select:
                # indicators, which span is going to be used
                self._axspan[trace] = ax[trace].axvspan(self._tmin,
                                                        self._tmax,
                                                        facecolor='y',
                                                        alpha=0.1)
                # selectors for timespan to be used
                self._spanselector[trace] = SpanSelector(ax[trace],
                                                         self._set_time_span,
                                                         'horizontal',
                                                         useblit=True,
                                                         minspan=self._minspan)

        # Add OK button
        # TODO: improve scaling of OK button size
        # self._ax_OK = ifigure.add_axes([ 0.9, 0.01, 0.06, 0.03])
        # button = plt.Button(self._ax_OK, 'OK')
        # def ap(event):
        #     self.close_ifig()
        # button.on_clicked(ap)

        ax[traces[-1]].set_xlabel('Time (s)')

        if description is None:
            ifigure.suptitle(self._pre_title)
        else:
            ifigure.suptitle(self._pre_title + ' for ' + description)

        return ifigure

    def update_ifig(self, timevector, data, samplingrate, traces, xlim=None):
        """
        Update the data of self.ifigure.

        Parameters
        ----------
        timevector : array
            The timeline of the data.
        data : 2D array
            The data to be plotted.
        samplingrate : float
            The samplingrate of the data in Hz. Data needs to have as many
            traces as the data that is currently plotted in `self.ifigure`.
        traces : str or list of str
            The traces contained in the data. `traces` has to contain the same
            traces as the data that is currently plotted in `self.ifigure`.
        xlim : tuple, optional
            Set the range (xmin, xmax) of the data to be plotted in units of
            timevector.
        """
        traces = tc.normalize(traces)
        ax = dict(list(zip(traces, self.ifigure.axes)))

        decimate = max(1, int(np.round(samplingrate / self.displayrate)))

        if xlim is None:
            xlim = (timevector[0], timevector[-1])

        # plot traces
        for i, trace in enumerate(traces):
            ax[trace].lines[0].set_data(timevector[::decimate],
                                        data[:, i][::decimate])
            ax[trace].set_xlim(*xlim)

        self.ifigure.canvas.draw()

    def close_ifig(self):
        """
        Close the figure the user has selected the region from and deregister
        all callback functions.
        """
        # remove references to SpanSelectors, to free memory
        self._axspan.clear()
        self._spanselector.clear()

        # remove onselect callbacks
        self._onselect_callbacks.clear()

        if self.ifigure is not None:
            # Remove button 'OK'
            # self._ax_OK.clear()
            # self._ax_OK.set_axis_off()
            # self._ax_OK = None

            # Some matplotlib backends will throw an error when trying to draw
            # the canvas. Simply ignoring the error that could happen here will
            # prevent the figure from not beeing closed, left open, and
            # preventing the next figure to be drawn.
            # Even though the except: pass clause is considered bad, here the
            # worst thing that could happen, is that the static figure produced
            # by the matplotlib backend upon closing is not up to date.
            # Therefore, except: pass should be considered as an acceptable
            # workaround for this case.
            try:
                # force redraw of the figure to remove the button in the nbagg
                # backend for the static figure
                self.ifigure.canvas.draw()

                # give backend time to redraw the figure (hack) before closing
                time.sleep(1)
            except:
                pass

            # close the figure
            plt.close(self.ifigure)

        self.ifigure = None

    def init_ifig(self, timevector, data, samplingrate, traces, tmin=None,
                  tmax=None, onselect_cb=None, xlim=None, description=None,
                  select=True):
        """
        Create and show a figure based on the given timevector and data that
        lets the user choose a timespan (tmin, tmax). Whenever the user changes
        the timespan, the function onselect_cb is called with parameters tmin
        and tmax.

        Parameters
        ----------
        timevector : np.array
            1D array with the timevector to plot
        data : np.array
            The data array containing the traces to plot
        samplingrate : float
            The samplingrate of the data in Hz
        traces : list of str
            List with names of traces to use for data, COLOR and NAME
        tmin : float, optional
            tmin for the span selector
        tmax : float, optional
            tmax for the span selector
        onselect_cb : function, optional
            A function which is called whenever the timespan was changed by the
            user. onselect_cb() is called with the two parameters tmin and
            tmax: onselect_cb(tmin, tmax).
        xlim : tuple, optional
            tuple with (tmin, tmax) to set the range of the data to plot.
            Defaults to (timevector[0], timevector[-1]).
        description : str, optional
            str describing the purpose of the region to select. Is used as a
            heading for the figure.
        """
        self.select = select

        # nbagg backend needs to have the figure closed and recreated
        # whenever the code of the cell displaying the figure is executed.
        # A simple update of the figure would let it disappear. Even a
        # self.ifigure.show() wouldn't work anymore.
        # For backends this just means a bit of extra calculation.
        # Therefore, close the figure first before replotting it.
        self.close_ifig()

        # replace old onselect_callback with new one
        self.set_onselect_callback(onselect_cb)

        # set plot parameters
        traces = tc.normalize(traces)
        plot_params = self._calculate_plot_params(len(traces), 1)
        set_plot_params(plot_params)

        # create interactive figure
        self.ifigure = self._ifigure(timevector,
                                     data,
                                     samplingrate,
                                     traces,
                                     tmin=tmin,
                                     tmax=tmax,
                                     xlim=xlim,
                                     description=description)

        # show the interactive figure
        self.ifigure.show()

    def _calculate_plot_params(self, rows, columns):
        """
        Calculate proper figure geometry and return matplotlib parameters
        accordingly. The geometry is calculated depending on the number of rows
        and columns.

        Parameters
        ----------
        rows : int
            Number of rows to plot data.
        columns : int
            Number of columns to plot data.

        Returns
        -------
        dict
            Matplotlib parameters. Can by set with `self.set_plot_params()`.
        """
        dpi = self._dpi

        fig_pixel_width  = columns * float(self._width)
        fig_pixel_height = rows    * float(self._height)

        left   =        self._left   / fig_pixel_width
        right  = 1.00 - self._right  / fig_pixel_width
        bottom =        self._bottom / fig_pixel_height
        top    = 1.00 - self._top    / fig_pixel_height

        wspace = self._wspace / float(self._width)
        hspace = self._hspace / float(self._height)

        plot_params = {
            'dpi': dpi,
            'fig_pixel_width': fig_pixel_width,
            'fig_pixel_height': fig_pixel_height,
            'left': left,
            'right': right,
            'bottom': bottom,
            'top': top,
            'wspace': wspace,
            'hspace': hspace
        }

        return plot_params

    def _set_time_span(self, tmin, tmax):
        """
        Set the timespan according to the SpanSelector and update the figure
        accordingly. Call all registered callback functions.

        This function is called upon any change of the SpanSelector.

        Parameters
        ----------
        tmin : float
            Start time point in s.
        tmax : float
            End time point in s.
        """
        self._tmin = tmin
        self._tmax = tmax
        for ax in list(self._axspan.values()):
                    ax.set_xy([[tmin, 0],  # lower left corner
                               [tmin, 1],  # upper left corner
                               [tmax, 1],  # upper right corner
                               [tmax, 0],  # lower right corner
                               [tmin, 0]])  # lower left corner
        self.ifigure.canvas.draw()
        self._fire_onselect_callbacks()

    def _fire_onselect_callbacks(self):
        """
        Call all registered callback functions.
        """
        # call onselect callback functions
        for function in self._onselect_callbacks:
            function(self._tmin, self._tmax)


def set_plot_params(plot_params=None):
    """
    Set matplotlib parameters.

    Parameters
    ----------
    plot_params : dict, optional
        Dictionary of matplotlib parameters to be set. The following default
        parameters are used, if not existent in plot_params:
        {
            'dpi' : 100.0,
            'fig_pixel_width' : 800.0,
            'fig_pixel_height' : 600.0,
            'left' : 65 / fig_pixel_width,
            'right' : 1.00 - 20 / fig_pixel_width,
            'bottom' : 55 / fig_pixel_height,
            'top' : 1.00 - 55 / fig_pixel_height,
            'wspace' : 18 / float(fig_pixel_width),
            'hspace' : 18 / float(fig_pixel_height),
            'lines.linewidth' : 0.25,
            'figure.facecolor' : '1.00',
            'figure.edgecolor' : '1.00',
            'axes.linewidth' : 1,
            'axes.labelsize' : 9,
            'legend.fontsize' : 9,
            'ytick.labelsize' : 7,
            'xtick.labelsize' : 7
        }
    """
    if plot_params is None:
        plot_params = {}

    dpi = plot_params.pop('dpi', 100.0)

    fig_pixel_width = plot_params.pop('fig_pixel_width', 800.0)
    fig_pixel_height = plot_params.pop('fig_pixel_height', 600.0)

    fig_width  = fig_pixel_width  / dpi
    fig_height = fig_pixel_height / dpi

    left   = plot_params.pop('left',         65 / fig_pixel_width)
    right  = plot_params.pop('right', 1.00 - 20 / fig_pixel_width)
    bottom = plot_params.pop('bottom',       55 / fig_pixel_height)
    top    = plot_params.pop('top',   1.00 - 55 / fig_pixel_height)

    wspace = plot_params.pop('wspace', 18 / float(fig_pixel_width))
    hspace = plot_params.pop('hspace', 18 / float(fig_pixel_height))

    lines_linewidth = plot_params.pop('lines.linewidth', 0.25)

    figure_facecolor = plot_params.pop('figure.facecolor', '1.00')
    figure_edgecolor = plot_params.pop('figure.edgecolor', '1.00')

    axes_linewidth = plot_params.pop('axes.linewidth', 1)

    axes_labelsize = plot_params.pop('axes.labelsize', 9)
    legend_fontsize = plot_params.pop('legend.fontsize', 9)
    ytick_labelsize = plot_params.pop('ytick.labelsize', 7)
    xtick_labelsize = plot_params.pop('xtick.labelsize', 7)

    rcparams = plt.rcParams.copy()
    rcparams = {
        'figure.figsize': (fig_width, fig_height),
        'figure.dpi': dpi, 'savefig.dpi': dpi,
        'figure.subplot.left': left,
        'figure.subplot.right': right,
        'figure.subplot.bottom': bottom,
        'figure.subplot.top': top,
        'figure.subplot.wspace': wspace,
        'figure.subplot.hspace': hspace,
        'figure.facecolor': figure_facecolor,
        'figure.edgecolor': figure_edgecolor,
        'lines.linewidth': lines_linewidth,
        'axes.linewidth': axes_linewidth,
        'axes.labelsize': axes_labelsize,
        'legend.fontsize': legend_fontsize,
        'ytick.labelsize': ytick_labelsize,
        'xtick.labelsize': xtick_labelsize
        }

    plt.rcParams.update(rcparams)


def browse_images(directory, prefix=None, suffix=None, extension='.png',
                  sort_key=None):
    """
    Display images of a directory with an interactive widget and a slider in a
    jupyter notebook, in the order sorted to their filename or a given key
    function.

    Parameters
    ----------
    directory : str
        The directory the images to be displayed are located in.
    prefix : str
        Display only the files beginning with `prefix`.
    suffix : str
        Display only the files ending with `suffix`.
    extension : str, optional
        The extension of the images that should be displayed. Default is
        '.png'.
    sort_key : function
        Function to be applied to every image filename found, before sorting.
    """
    images = hp.files(directory, prefix, suffix, extension, sort_key)
    stop = len(images)
    if stop < 1:
        print("No file found with prefix '%s', suffix '%s', and extension '%s'"
              % (prefix, suffix, extension))
        return

    def view_image(i):
        image = images[i]
        try:
            display(Image(filename=image))
            print(image)
        except:
            print('No image found!')

    slider = IntSlider(min=0, max=stop-1, step=1, value=0,
                       description='Image:')

    return interact(view_image, i=slider)
