# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:22:31 2016

@author: Tobias Jachowski
"""
import collections
import matplotlib.pyplot as plt
import numpy as np
from abc import ABCMeta, abstractmethod

from .. import gui
from .. import helpers as hp
from .. import traces as tc
from ..evaluate import signal as sn
from ..graph import GraphMember
from ..picklable import InteractiveAttributes


class GraphicalMod(object):
    """
    This class's subclasses should implement `_figure()` and `_update_fig()`,
    which return and update a matplotlib figure, respectively. The figure can
    be accessed by `self.figure`.

    Parameters
    ----------
    figure
    modification : Modification
    """
    def __init__(self, modification=None, **kwargs):
        # Register the modification which should be graphically adjusted
        self.modification = modification

        # Initialize figure to None, which effectively disables
        # `self.update_fig()` and Co. and prevent them from throwing an error
        self._fig = None

    def _set_plot_params(self, plot_params=None):
        if plot_params is None:
            plot_params = {}
        gui.set_plot_params(plot_params=plot_params)

    def display(self, plot_params=None):
        self.init_fig(plot_params=plot_params)

    def init_fig(self, show=True, plot_params=None):
        """
        This method calls self._figure() to create an interactive figure and
        interact with the user to determine the parameters necessary to
        calculate the modification (see self._recalculate()). and
        self._close_fig() to release all references to the actors of the
        figure.
        `self._figure()` and self._close_fig() should be (over)written by
        subclasses.
        """
        # Only create a figure, if the function `self._figure()` is implemented
        if not hasattr(self, '_figure'):
            return

        # close the figure
        # nbagg backend needs to have the figure closed and recreated
        # whenever the code of the cell displaying the figure is executed.
        # A simple update of the figure would let it disappear. Even a
        # self.figure.show() wouldn't work anymore.
        # For backends this just means a bit of extra calculation.
        # Therefore, close the figure first before replotting it.
        self.close_fig()

        # set default plot parameters, can be recalled / overwritten in
        # `self._figure()`
        self._set_plot_params(plot_params=plot_params)

        # create the figure
        self.figure = self._figure()

        # update the figure
        self.update_fig()

        # show the figure
        if show:
            self.figure.show()

    def update(self, **kwargs):
        self.update_fig(**kwargs)

    def update_fig(self, **kwargs):
        if self._fig is not None:
            self._update_fig(**kwargs)
            self._figure_canvas_draw()

    def _update_fig(self, **kwargs):
        pass

    def close_fig(self):
        if self._fig is not None:
            self._pre_close_fig()
            self._close_fig()
            self._post_close_fig()

    def _pre_close_fig(self):
        """
        Method to be overwritten by subclasses.
        """
        pass

    def _close_fig(self):
        # force redraw of the figure
        self._figure_canvas_draw()

        # close the figure
        plt.close(self.figure)

        # release memory
        self.figure = None

    def _post_close_fig(self):
        """
        Method to be overwritten by subclasses.
        """
        pass

    def _figure_canvas_draw(self):
        # Some matplotlib backends will throw an error when trying to draw the
        # canvas. Simply ignoring the error that could happen here will prevent
        # the figure from not beeing closed, left open, and preventing the next
        # figure to be drawn. Even though the "except: pass" clause is
        # considered bad, here the worst thing that could happen is that the
        # figure produced by the matplotlib backend upon closing is not
        # updated. Therefore, "except: pass" should be considered as an
        # acceptable workaround for this case.
        try:
            # redraw the figure, before closing it
            self.figure.canvas.draw()
        except:
            pass

    @property
    def figure(self):
        """
        The matplotlib figure that represents and/or adjusts the parameters of
        `self.modification`.
        """
        # Automatically initialize a figure
        if self._fig is None:
            self.init_fig(show=False)
        # Return a previously initialized figure
        return self._fig

    @figure.setter
    def figure(self, figure):
        self._fig = figure


class Modification(GraphMember, metaclass=ABCMeta):
    """
    Modification is an abstract class, that implements methods to modify the
    data of a `View` (`view_apply`) and adjust the parameters which control the
    behaviour of the modifications applied.

    Whenever one of the parameters needed to calculate the modification is
    changed, the view, this modification is applied to, is informed.
    `self.set_changed()` Has to be called upon any change of the modification
    that influences the behaviour of `self.modify()`. In essence, these are all
    parameters that are used to determine the modification. Therefore, this
    should be called by all setters of the parameters/attributes.

    Every subclass of Modification has to implement a constructor method
    `self.__init__(self, **kwargs)`, which calls the superclasses' constructor
    and sets the traces, the modification is applied to with the keyword
    parameter `traces_apply`. An example could be:
    super().__init__(traces_apply=['psdX', 'psdZ'], **kwargs)
    """
    # set a graphical modification, which will, per default, do nothing
    GRAPHICALMOD = GraphicalMod

    def __init__(self, traces_apply=None, view_apply=None, view_based=None,
                 automatic_switch=False, datapoints=-1, **kwargs):
        # Call the constructor of the superclass `GraphMember` and set the
        # maximum allowed number of parents (`view_based`) and childs
        # (`view_apply`) to one.
        super().__init__(max_children=1, max_parents=1, **kwargs)

        # A `Modification` has to be applied to a `View`!
        if view_apply is None:
            raise TypeError("Modification missing required positional argument"
                            " `view_apply`.")

        # Set the view, from where the parameters for the modification are
        # calculated from
        if view_based is not None:
            self.view_based = view_based

        # Set the view, whose data is going to be modified
        self.view_apply = view_apply

        # Set the traces, which are modified by this `Modification`
        self.traces_apply = traces_apply

        # Initialize InteractiveAttributes object, which will hold all the
        # parameters that the user should interact with.
        self.iattributes = InteractiveAttributes()

        # A checkbox to switch on/off the automatic determination of the
        # parameters that are used to calculate the modification in the method
        # `self.recalculate()`. The attribute `self.automatic` is checked in
        # the method `self.recalculate()`. If `automatic` is True, the
        # parameters are recalculated, otherwise the parameters are left
        # unchanged. Whenever `automatic` is changed (by the user or
        # automatically), `self.evaluate()` is called.
        if automatic_switch:
            self.add_iattribute('automatic', description='Automatic mode',
                                value=True, unset_automatic=False,
                                set_changed=False,
                                callback_functions=[self.evaluate])

        # A checkbox to de-/activate this `Modification`. This attribute gets
        # evaluated by `self.modify()`. If the `Modification` is active, it
        # modifies data, otherwise not, i.e. modify() returns modified or
        # unmodified original data, respectively.
        desc = "".join((self.__class__.__name__, " active"))
        self.add_iattribute('active', description=desc, value=True,
                            unset_automatic=False)

        # Datapoints is used to calculate and/or present modification. The
        # attribute `datapoints` is used to calculate a decimating factor and
        # speed up the calculations and/or plot commands.
        if datapoints > 0:
            desc = "Datapoints to calculate/visualize modification"
            self.add_iattribute('datapoints', description=desc,
                                value=datapoints, unset_automatic=False)

        # Add a Button to manually call the method `self.evaluate()`.
        self.add_iattribute('evaluate', description='Evaluate',
                            unset_automatic=False, set_changed=False,
                            callback_functions=[self.evaluate])

    def add_iattribute(self, key, description=None, value=None,
                       unset_automatic=True, set_changed=True,
                       callback_functions=None, **kwargs):
        """
        Add logic for automatic checkbox.
        Register widget with unset_automatic=True
        (-> Upon change of widget, unset automatic mode).
        Change default behaviour by setting kwarg: unset_automatic = False
        Add logic for triggering changed (calling self.set_changed).
        Register widget with set_changed=True.
        """
        if callback_functions is None:
            callback_functions = []

        if unset_automatic:
            callback_functions.append(self._unset_automatic)

        if set_changed:
            callback_functions.append(self.set_changed)

        self.iattributes.add(key, description=description, value=value,
                             callback_functions=callback_functions, **kwargs)

    def _unset_automatic(self, leave_automatic=False, **kwargs):
        """
        Add the logic for the automatic checkbox. If the value of an attribute
        is changed and the attribute was created with `unset_automatic=True`,
        deactivate the automatic mode (see `self.add_iattribute()`). To
        temporarily leave the automatic mode status untouched when changing the
        value of an attribute, i.e. not unset the automatic mode, set the value
        of the attribute with the keyword argument `leave_automatic=True`
        (see method `self.iattributes.set_value()`)
        """
        if not leave_automatic:
            self.iattributes.set_value('automatic', False, callback=False)

    def evaluate(self):
        """
        Implement the (re)calculation for the values necessary to calculate the
        modification in the subclass and call recalculate() of the superclass
        (this class).
        """
        if self.updated:
            # This method makes sure the modification is calculated with the
            # current values of the View this modification is based on. It is
            # called by self.modify().

            # When a View requests data, it calls modify(), which in turn calls
            # recalculate(). Recalculate(), if necessary, calls
            # get_data_modified() from the View it is based on, which again
            # triggers a call of modify() and a subsequent recalcaulte() of all
            # modifications associated with this View.
            # Modification need update, because view, this mod is based on,
            # was changed.
            # self._view_based.evaluate()is not needed, it is called via:
            # recalculate() -> get_data_based() -> _view_based.get_data() ->
            # get_modified_data() -> super().evaluate()
            return

        # Recalculate and print info of recalculated values if in automatic
        # mode
        if self.recalculate():
            self.print_info()

        # Update figure after recalculation has taken place
        self.graphicalmod.update()

    def recalculate(self):
        # Check if recalculation of parameters is necessary
        if self.updated:
            return False
        # Check the attribute self.automatic, whether the parameters needed for
        # the calculation of the modification should be determined
        # automatically or not. If values are set manually, no recalculation is
        # necessary, and `self` is therefore up to date.
        if not self.automatic:
            self.updated = True
            return True
        # Recalculate the parameters, inform the view this `Modification`
        # is applied to about the change, and set `self` to be updated.
        self._recalculate()
        self.set_changed(updated=True)
        return True

    def _recalculate(self):
        """
        This method should be overwritten by subclasses and perform the
        recalculation necessary to determine the parameters used by this
        Modification to modify the data in `self._modify()`.
        """
        pass

    def print_info(self):
            print("Values for Modification of class %s:"
                  % self.__class__.__name__)
            if not self.automatic:
                print("  Parameters set manually!")
            for key, widget in self.iattributes._widgets.items():
                if hasattr(widget, 'value'):
                    if isinstance(widget.value, float):
                        print("    %s: %.5f" % (widget.description, widget.value))
                    if isinstance(widget.value, collections.Iterable):
                        print("    %s: %s" % (widget.description, widget.value))
            self._print_info()

    def _print_info(self):
        """
        This method should be overwritten by subclasses, which want to print
        extra info additionally to the info of the calculated paremeters.
        """
        pass

    def modify(self, data, samples, traces_idx):
        """
        Modifies data and returns the modified array.

        Parameters
        ----------
        data : 2D numpy.ndarray of type float
            `data` holds the data to be modified
        samples : index array or slice
            `samples` is the index of the samples that was used to get the
            `data`
        traces : index array or slice
            `traces` is the index of the traces that was used to get the `data`
        """
        # Modification is active.
        if self.active:
            # Check if traces contained in data are modified by this
            # modification.
            data_traces = self.view_apply.idx_to_traces(traces_idx)
            mod_traces = self.traces_apply

            # Calculate the indices of traces contained in data and
            # modification. First, calculate indices of modification traces.
            mod_index = hp.overlap_index(mod_traces, data_traces)
            if len(mod_index) > 0:
                # At least one trace exists in both data and modification.
                # Therefore, the data needs to be modified...
                mod_index = hp.slicify(mod_index)

                # Calculate indices of traces of the data in such a way that
                # `data[:, data_index]` indexes the same traces as
                # `self.traces_apply[mod_index]`
                data_index = np.array([data_traces.index(trace)
                                       for trace
                                       in np.array(mod_traces)[mod_index]])
                data_index = hp.slicify(data_index)

                # Trigger a recalculation of the parameters for the
                # modification (if necessary) before modifying the data.
                self.evaluate()

                # Modify and return the modified data
                return self._modify(data=data,
                                    samples=samples,
                                    data_traces=data_traces,
                                    data_index=data_index,
                                    mod_index=mod_index)
        # Return unmodified data
        return data

    @abstractmethod
    def _modify(self, data, samples, data_traces, data_index, mod_index):
        """
        Is called by self.modify() whenever data is requested and needs to be
        modified.

        Parameters
        ----------
        data : 2D numpy.array()
            Contains the data, indexed by samples and data_traces
        samples : slice or 1D numpy.array()
            Is the index of the samples contained in data, which was
            given/asked by the user/process who called _get_data().
        data_traces : list of str
            Contains a list of traces (str) existent in data, which
            was given/asked by the user/process who called _get_data().
        data_index : slice or 1D numpy.array()
            data[:, data_index] gives the data, which is modified by
            this modification
        mod_index : slice or 1D numpy.array()
            np.array(self.traces_apply)[mod_index] gives the traces,
            which are existent in data and also modified by this modfication.

        Returns
        -------
        2D numpy.array()
            The modified data.
        """
        # modify data here, like so:
        # data[:,data_index] -= modification[:,mod_index]
        return data

    @property
    def updated(self):
        return self._updated

    @updated.setter
    def updated(self, value):
        """
        Gets set to True, after all `Views`, this `Modification` is based on,
        have been updated and after this `Modification` has been recalculated.
        This is automatically taken care of by `self.evaluate()` ->
        `self.recalculate()`.

        Gets called by a `View`, this `Modification` is based on, whenever the
        `View` (a `Modification` of the `View`) has been changed. It
        automatically informs its own `View`, that there was a change, by
        calling `self.set_changed()`.
        """
        self._updated = value

    def member_changed(self, ancestor=True, calledfromself=False,
                       index_shift=None, **kwargs):
        # If a change of an ancestor View or a MultiRegion was triggered by an
        # index_shift, the modification needs to recalculate itself, i.e.
        # the modification will alter its changeing behaviour. Because an
        # index_shift change is only transmitted to `level=1`, inform the
        # descendants of the change itself. A change of descendants is ignored.
        if index_shift is not None and not calledfromself and ancestor:
            self.set_changed(includeself=False)

        # Update update status
        super().member_changed(ancestor=ancestor,
                               calledfromself=calledfromself, **kwargs)

    def _get_data(self, based=True, samples=None, traces=None, window=False,
                  decimate=False, copy=True):
        if based:
            view = self.view_based
        else:
            view = self.view_apply

        if not isinstance(window, bool) and isinstance(window, int):
            window = window
        elif window:
            window = self.decimate
        else:
            window = 1

        if not isinstance(decimate, bool) and isinstance(decimate, int):
            decimate = decimate
        elif decimate:
            decimate = self.decimate
        else:
            decimate = 1

        if not based:
            old_active = self.iattributes.active
            self.iattributes.set_value('active', False, callback=False)

        data = view.get_data(traces=traces, samples=samples,
                             moving_filter='mean', window=window,
                             decimate=decimate, copy=copy)

        if not based:
            self.iattributes.set_value('active', old_active, callback=False)

        return data

    def _get_data_based(self, samples=None, traces=None, window=False,
                        decimate=False, copy=True):
        """
        decimate is False per default. If decimate is True, it only gets used,
        if samples are set to None (step information in samples precedes over
        decimate).
        """
        return self._get_data(based=True, samples=samples, traces=traces,
                              window=window, decimate=decimate, copy=copy)

    def _get_data_apply(self, samples=None, traces=None, window=False,
                        decimate=False, copy=True):
        """
        Get data of view apply with all modifications applied, except self.
        This is achieved by setting the self.__active flag to False.
        self.__active is intentionally set directly by accessing the attribute
        and not using the property/set_active() method, to prevent firing the
        self.set_changed() method within the set_active() method.
        decimate is False per default. If decimate is True, it only gets used,
        if samples are set to None (step information in samples precedes over
        decimate).
        """
        return self._get_data(based=False, samples=samples, traces=traces,
                              window=window, decimate=decimate, copy=copy)

    def calculate_bin_means(self, data=None, traces=None, binnumber=None,
                            datapoints_per_bin=None, sorttrace=0):
        """
        Calculates binned means based on the data to be fitted. The binned
        means are usually used by data fitting routines.

        Parameters
        ----------
        data : 2D numpy.ndarray of type float, optional
            Defaults to `self._get_data_based(traces=traces, decimate=True)`.
        traces : str or list of str, optional
            Defaults to `self.traces_apply`.
        binnumber : int, optional
            Number of bins that contain the datapoints to be averaged. If
            possible, it defaults to (`self.iattributes.datapoints` /
            `datapoints_per_bin`), otherwise binnumber defaults to
            (`self.view_based.datapoints` / `datapoints_per_bin`).
        datapoints_per_bin : int, optional
            Average number of datapoints to be averaged in one bin. Defaults to
            25.
        sorttrace : int, optional
            Trace (column) of `data` that acts as sorting index upon binning
            for the rest of the data. Defaults to the first trace (column) in
            the data.

        Returns
        -------
        1D numpy.ndarray of type float
            The averaged bin values.
        float
            The size of one bin.
        """
        # Bin data and average bins to prevent arbitrary weighting of bins with
        # more datapoints
        if binnumber is None:
            binnumber = self._binnumber(datapoints_per_bin=datapoints_per_bin)

        # get the traces to retrieve data from
        if traces is None:
            traces = self.traces_apply

        # get the data to bin
        if data is None:
            data = self._get_data_based(traces=traces, decimate=True)

        # create the bins based on one trace of the data
        minimum = np.min(data[:, sorttrace])
        maximum = np.max(data[:, sorttrace])
        bins = np.linspace(minimum, maximum, binnumber + 1)

        # Get the indices of the bins to which each value in input array
        # belongs.
        indices = np.digitize(data[:, sorttrace], bins)

        # fill the bins with the means of the data contained in each bin
        bin_means = np.array([data[indices == i].mean(axis=0)
                              for i in range(1, len(bins))
                              if np.any(indices == i)])

        bin_size = bins[1] - bins[0]

        return bin_means, bin_size

    def _binnumber(self, datapoints_per_bin=None):
        # On average 25 datapoints per bin
        datapoints_per_bin = datapoints_per_bin or 25
        if 'datapoints' in self.iattributes:
            binnumber = self.iattributes.datapoints / datapoints_per_bin
        else:
            binnumber = self.view_based.datapoints / datapoints_per_bin
        binnumber = max(1, int(np.round(binnumber)))
        return binnumber

    _NAME = {
        'position': ['positionX', 'positionY'],
        'psd': ['psdX', 'psdY'],
        'axis': ['X', 'Y']
    }

    def _excited(self, traces=None):
        traces = traces or ['positionX', 'positionY']
        data = self._get_data_based(traces=traces, copy=False)
        return sn.get_excited_signal(data)

    def interact(self):
        self.recalculate()
        self.iattributes.display()
        self.graphicalmod.display()

    @property
    def graphicalmod(self):
        # ZODB volatile
        if not hasattr(self, '_v_graphicalmod'):
            self._v_graphicalmod \
                = self.__class__.GRAPHICALMOD(modification=self)
        return self._v_graphicalmod

    @property
    def active(self):
        active = False
        if 'active' in self.iattributes:
            active = self.iattributes.active
        return active

    @property
    def automatic(self):
        # Does the modification automatically calculate its parameters
        automatic = True
        if 'automatic' in self.iattributes:
            automatic = self.iattributes.automatic
        return automatic

    @property
    def datapoints(self):
        if 'datapoints' in self.iattributes:
            return self.iattributes.datapoints
        else:
            return self.view_based.datapoints

    @property
    def decimate(self):
        if 'datapoints' in self.iattributes:
            return max(1, int(np.round(self.view_based.datapoints
                                       / self.datapoints)))
        else:
            return 1

    @property
    def view_based(self):
        return self.parent

    @property
    def view_apply(self):
        return self.child

    @view_based.setter
    def view_based(self, view):
        self.set_parent(view)

    @view_apply.setter
    def view_apply(self, view):
        self.set_child(view)

    def lia(self, trace):
        """
        Return the local index of trace in traces_apply
        """
        return self.traces_apply.index(trace)

    @property
    def traces_apply(self):
        # return a copy to protect local copy
        return self._traces_apply.copy()

    @traces_apply.setter
    def traces_apply(self, traces):
        if traces is None:
            traces_apply = []
        else:
            traces_apply = tc.normalize(traces)
        self._traces_apply = traces_apply
