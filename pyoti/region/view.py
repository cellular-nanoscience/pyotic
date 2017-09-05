# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 22:05:56 2016

@author: Tobias Jachowski
"""
import collections
import numpy as np
import operator

from .region import Region
from .record import Record
from ..modification import Modification


class View(Region):
    """
    Defines a View on a given Region (Record or View), the parent of that view.
    The parent is an instance of a MultiRegion, which gives the options to add
    even further parents (see MultiRegion).
    The region of a View of the data to refer to can (in contrast to a Record)
    be choosen by setter methods (start, stop, tmin, tmax).
    """

    def __init__(self, parent=None, sticky_start=True, sticky_stop=True,
                 **kwargs):
        """
        Parameters
        ----------
        parent : Region
        sticky_start : bool
        sticky_stop : bool
        **kwargs

        Attributes
        ----------
        records
        record
        views
        modifications_apply
        modifications_based
        group
        parent
        calibration
        num_traces
        traces
        samplingrate
        start
        stop
        sticky_start : bool
            If True and `self.start` has the same position as the parent's
            start, `self.start` follows parent's start automatically, if it
            is changed.
            If False, only shrinkage of the parent's region (increasing start)
            updates `self.start`, if necessary.
            `self.sticky_start` is automatically set to False, if `self.start`
            is set manually.
        sticky_stop : bool
            If True and `self.stop` has the same position as the parent's
            stop, `self.stop` follows parent's stop automatically, if it
            is changed.
            If False, only shrinkage of the parent's region (decreasing stop)
            updates `self.stop`, if necessary.
            `self.sticky_stop` is automatically set to False, if `self.stop`
            is set manually.
        """
        super().__init__(**kwargs)

        if parent is None:
            raise TypeError("View missing 1 required positional argument: "
                            "`parent`")
        # Create a new MultiRegion with `parent` as parent and add it to `self`
        # as parent
        multiregion = MultiRegion(parent=parent, name=self.name)
        if not self.add_parent(multiregion):
            raise ValueError("Could not add `self` as a child, most probably "
                             "circular reference detected.")

        self.start = 0
        self.stop = self.parent.datapoints

        # Determine, whether start and stop (tmin and tmax) should follow the
        # increase of the parent region or compensate for it to point to the
        # same region.
        self.sticky_start = sticky_start
        self.sticky_stop = sticky_stop

    @property
    def records(self):
        return self.parent.records

    @property
    def record(self):
        return self.parent.record

    @property
    def views(self):
        return self.parent.views

    @property
    def modifications_apply(self):
        return self.modifications(key='apply')

    @property
    def modifications_based(self):
        return self.modifications(key='based')

    def modifications(self, key='apply'):
        if key == 'apply':
            return self.parent_instances(Modification)
        else:
            return self.child_instances(Modification)

    @property
    def group(self):
        modifications = list(self.modifications_based)
        if len(modifications) >= 1:
            # one or more modification_based -> same group as first
            # modification
            return modifications[0].group
        else:
            # no modification_based -> own group, which could be the same but
            # doesn't have to
            return self._group

    @group.setter
    def group(self, group):
        modifications = list(self.modifications_based)
        if len(modifications) >= 1:
            # one or more modification_based -> set group of first modification
            modifications[0].group = group
        else:
            # no modification_based -> set own group
            # self.group = None (or given group) is called in
            # super().__init__(), which ensures initialisation of self._group
            self._group = group

    @property
    def parent(self):
        # MultiRegion is the parent of View
        return self.first_ancestor_instance(MultiRegion, level=1, dft=False)

    @property
    def calibration(self):
        return self.parent.calibration

    @property
    def num_traces(self):
        if self.parent is None:
            return 0
        return self.parent.num_traces

    @property
    def traces(self):
        if self.parent is None:
            return []
        return self.parent.traces

    @property
    def samplingrate(self):
        if self.parent is None:
            return 0.0
        return self.parent.samplingrate

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @start.setter
    def start(self, start):
        if not hasattr(self, '_start'):
            # Initialise start, when called for the first time
            self._start = start
        else:
            if start != self._start:
                self.sticky_start = False
                # Validate start to be within 0 and self.parent.datapoints
                start = max(start, 0)
                start = min(start, self.parent.datapoints)

                # create shift for auto shift update of children
                shift = (self, 'start', 0, start - self._start)

                self._start = start

                # inform children - views and modifications based on this
                # view about change
                self.set_changed(level=1, index_shift=shift)

    @stop.setter
    def stop(self, stop):
        if not hasattr(self, '_stop'):
            # Initialise stop, when called for the first time
            self._stop = stop
        else:
            if stop != self._stop:
                self.sticky_stop = False
                # Validate stop to be within 0 and self.parent.datapoints
                stop = max(stop, 0)
                stop = min(stop, self.parent.datapoints)

                # create shift for auto shift update of children
                shift = (self,
                         'stop',
                         self._stop - self._start,
                         stop - self._stop)

                self._stop = stop

                # inform children - views and modifications based on this
                # view about change
                self.set_changed(level=1, index_shift=shift)

    def member_changed(self, ancestor=True, calledfromself=False,
                       index_shift=None, **kwargs):
        # If a change of an ancestor View or a MultiRegion was triggered by an
        # index_shift, check, if this index_shift of the ancestor necessitates
        # an update of self.stop or self.start. If yes, then set self.stop or
        # self.start accordingly and implicitly inform further descendants of
        # the change (The children of a View or a MultiRegion are informed upon
        # an index_shift, see call of set_changed() in property start/stop).
        # A change of descendants should be ignored.
        if index_shift is not None and not calledfromself and ancestor:
            if not self._update_index(index_shift):
                # It was not necessary to change start and stop, this means,
                # the data this View offers was not changed and in turn the
                # cached data does not have to be recalculated -> return and do
                # not trigger the recaching nor set the self.updated to False
                # (-> do not call super.member_changed())
                return

        # Set or unset cached data and update update status
        super().member_changed(ancestor=ancestor,
                               calledfromself=calledfromself, **kwargs)

    def _update_index(self, index_shift):
        """
        If the parent has changed its indexes, check, whether the own indexes
        need to be adjusted, too and inform child about index_shift, if
        necessary.

        Parameters
        ----------
        index_shift : tuple of (Region, str, int, int)
            The tuple contains: The caller Region that caused the index shift,
            the name of which index was changed (start/stop), the original
            index of the position where the change took place and the amount by
            which start or stop were shifted (caller, name, index, shift).

        Returns
        -------
        bool
            When `self.start` or `self.stop` had to be updated, return True,
            otherwise False.
        """
        # Get the start/stop shift from the parent
        caller, name, index, shift = index_shift
        # print("%s -> %s" %(caller.name, self.name))
        # print("%s, %s, %s" % (name, index, shift))

        # Initialize start and stop, in case there is now shift
        start = self._start
        stop = self._stop

        returnvalue = False

        if name == 'start':

            # There are three possible cases for self._start:
            # 1. start > index
            #    change start
            #    if shift > start - index:
            #      change start AND inform child
            # 2. start == index
            #    if shift is positive:
            #      leave start, inform child
            #    if shift is negative:
            #      if sticky: leave start, inform child
            #      else: change start
            # 3. start < index
            #    leave start, inform child
            #
            # There is one case, where self._stop is changed:
            # 1. stop > index
            #    change stop

            # Change self._start
            if self.sticky_start:
                op = operator.gt  # start > index is sticky
            else:
                op = operator.ge  # start >= index is non sticky
            if op(self._start, index):
                # Validate start to be at least as great as start of parent
                start = max(self._start - shift, index)

            # Change self._stop
            # and keep self._stop > index (start)
            if self._stop > index:
                # Validate stop to be at min as great as start + 1 of parent
                stop = max(self._stop - shift, index + 1)

            # Calculate index_shift for child
            _index = max(index - self._start, 0)
            _shift = max(shift - max(self._start - index, 0), 0)
            if self._start == index and self.sticky_start:
                _shift = shift

            # Set changed start and stop
            if self._start != start or self._stop != stop:
                # print(" Setting start in %s to %s" % (self.name, start))
                # print(" Setting stop in %s to %s" % (self.name, stop))
                # Set start and/or stop
                self._start = start
                self._stop = stop
                returnvalue = True

            # Inform child about 'start' index shift, if self._start and
            # self._stop do not compensate for the index shift
            if _shift != 0:
                i_shift = (self, name, _index, _shift)
                self.set_changed(level=1, index_shift=i_shift,
                                 includeself=False)

        else:  # name == 'stop':

            # Change self._start
            if self._start >= index:
                # This start follows the parent's stop index change
                start = self._start + shift
            else:  # self._start < self.parent.
                # Keep self._start < index (stop)
                # Validate start to be at max as great as stop - 1 of parent
                start = min(self._start, index + shift - 1)

            # Change self._stop
            if self.sticky_stop:
                op = operator.ge  # >= is sticky
            else:
                op = operator.gt  # > is non sticky
            if op(self._stop, index):
                # This start follows the parent's stop index change
                stop = self._stop + shift
            else:  # self._stop still could be after parent's start, therefore:
                stop = min(self._stop, index + shift)

            # Calculate index_shift for child
            _index = max(index - self._start, 0)
            _shift = min(shift + max(index - self._stop, 0), 0)
            if self._stop == index and self.sticky_stop:
                _shift = shift

            # Set changed start and stop
            if self._start != start or self._stop != stop:
                # print(" Setting start in %s to %s" % (self.name, start))
                # print(" Setting stop in %s to %s" % (self.name, stop))
                # Set start and/or stop
                self._start = start
                self._stop = stop
                returnvalue = True

            # Inform child about 'stop' index shift, if self._start and
            # self._stop do not (fully) compensate for the index shift
            if _index != 0 and _shift != 0:
                i_shift = (self, name, _index, _shift)
                self.set_changed(level=1, index_shift=i_shift,
                                 includeself=False)

        return returnvalue

    def _get_data_uncached(self, samples, traces_idx, copy):
        """
        Return uncached data, i.e.:
        Correct samples such that they relate to self.parent.datapoints, rather
        than self.datapoints, i.e. shift index by self.start.
        Get data from parent.
        Apply all modifications registered at this View.
        Return data.
        """

        # parentize the requested samples, i.e. correct for self.start
        if isinstance(samples, slice):
            p_start = samples.start + self.start
            p_stop = samples.stop + self.start
            p_step = samples.step
            p_samples = slice(p_start, p_stop, p_step)
        else:
            # samples is an np.array
            p_samples = samples + self.start

        data = self.parent._get_data(p_samples, traces_idx, copy)

        # modify data by applying all modifications to the unmodified data
        for mod in self.modifications():
            data = mod.modify(data, samples, traces_idx)

        return data


class MultiRegion(Region):
    """
    Holds references to multi parent instances (Region), offers methods
    to add and remove a parent and a method to serve the data of all parents
    concatenated. Inherits from class Region.
    """
    def __init__(self, parent=None, **kwargs):
        """
        Parameters
        ----------
        parent : Region or Iterable of Regions
        **kwargs

        Attributes
        ----------
        group
        start
        stop
        records
        views
        record
        calibration
        num_traces
        traces
        samplingrate
        """
        super().__init__(max_children=1, **kwargs)

        if parent is None:
            raise TypeError("Parent missing 1 required positional argument: "
                            "'parent'")

        if isinstance(parent, collections.Iterable):
            for p in parent:
                self.add_parent(p)
        else:
            self.add_parent(parent)

    def member_changed(self, ancestor=True, calledfromself=False,
                       index_shift=None, **kwargs):
        # Possible causes of an index shift:
        #    1. index shift of parent
        #    2. adding or removing a parent
        # This method is called in case 1. For handling of case 2 see
        # methods self.add_parent() and self.remove_parent().
        #
        # If a change of an ancestor View or a MultiRegion was triggered by an
        # index_shift, check, if this index_shift of the ancestor necessitates
        # an update of self.stop or self.start. If yes, then set self.stop or
        # self.start accordingly and implicitly inform further descendants of
        # the change (only the children of a View or a MultiRegion (level=1)
        # are informed upon an index_shift, see call of set_changed() in
        # property start/stop).
        # A change of descendants should be ignored.
        if index_shift is not None and not calledfromself and ancestor:
            caller, name, index, shift = index_shift
            # print("%s -> %s" %(caller.name, self.name))
            # print("%s, %s, %s" % (name, index, shift))

            # Find the start index of the changed caller
            start = 0
            for parent in self.parents:
                if parent is caller:
                    break
                start += parent.datapoints

            # Adjust index according to start index of parent/caller
            index = start + index

            # Create new index_shift based on old one
            index_shift = (self, name, index, shift)

            # Inform children - usually views - based on this MultiRegion about
            # index shift change
            self.set_changed(level=1, index_shift=index_shift,
                             includeself=False)

        # Set or unset cached data and update update status
        super().member_changed(ancestor=ancestor,
                               calledfromself=calledfromself, **kwargs)

    @property
    def group(self):
        try:
            return self.child.group
        except:
            return None

    @group.setter
    def group(self, group):
        try:
            self.child.group = group
        except:
            pass

    @property
    def start(self):
        return 0

    @property
    def stop(self):
        # the sum of all parent.datapoints
        datapoints = 0
        for parent in self.parents:
            datapoints += parent.datapoints
        return datapoints

    @property
    def records(self):
        """
        Return all parents of instance Record as a generator.
        """
        return self.parent_instances(Record)

    @property
    def views(self):
        """
        Return all parents of instance View as a generator.
        """
        return self.parent_instances(View)

    @property
    def record(self):
        """
        Returns the first record that view and subsequently that view's parents
        are based on.
        """
        # TODO: implement better decision, which record to return, depending on
        # the parents; this decision influences the properties: calibration,
        # samplingrate, and num_traces
        # Presently, only views/records can be added as parent, whose
        # calibration, samplingrate, and num_traces are equal to this view
        return self.first_ancestor_instance(Record)

    @property
    def calibration(self):
        # TODO: Views combine the raw signals -> no handling of different
        # calibration in View meaningful -> Create higher grouping system which
        # gives access to calibrated distances/extensions/forces, only, but not
        # the raw signal -> groups of Results? -> shift Calibration of Record
        # to Reults
        return self.record.calibration

    @property
    def num_traces(self):
        if self.record is None:
            return 0
        # TODO: Implemement handling of different num_traces: Allow only the
        # overlap of different traces to be used?
        return self.record.num_traces

    @property
    def traces(self):
        if self.record is None:
            return []
        return self.record.traces

    @property
    def samplingrate(self):
        if self.record is None:
            return 0.0
        # TODO: Implement handling of different samplingrate: up/down sampling?
        # Pandas? Every parent has its own samplingrate.
        return self.record.samplingrate

    def add_parent(self, region, index=None, after=None, before=None,
                   update_indices=True):
        """
        Priority of inserting the new region at a specific index in descending
        order:

            #. index
            #. after
            #. before

        Returns
        -------
        bool
            True if parent could be added or parent is already present in
            self._parents. Otherwise return False.
        """
        if region is None:
            return False

        # TODO: Implement possibility of handling different calibrations,
        # samplingrate, and num_traces
        # Check, if there is already at least one view existent and if so,
        # make sure the new parent is based on the same calibration, has
        # the same samplingrate and the same number of traces
        if (next(self.parents, None) is not None
            and (self.calibration is not region.calibration
                 or self.samplingrate != region.samplingrate
                 or self.num_traces != region.num_traces)):
            return False

        # Possible causes of an index shift:
        #    1. index shift of parent
        #    2. adding or removing a parent
        # This method is called in case 2. For handling of case 1 see method
        # self.member_changed().

        # Add parent.
        # Adding a parent usually causes a change of the children. Therefore,
        # super().add_parent() would call set_changed() itself, without any
        # parameters. This means, if a parent was added, set_changed() would
        # trigger a recalculation of Modification, regardless of whether the
        # change affected the region selected by the child View, or not.
        # However, adding a parent to a View only causes an index shift.
        # Therefore, disable generic set_changed() mechanism and initiate a
        # specific one with the index_shift information.
        if region in self.parents:
            # Region is already in parents.
            # Avoid informing the children about index_shift.
            # This is necessary, because further down super().add_parent()
            # would return True, and would therefore not prevent the
            # information about the index_shift.
            print("Region %s already referenced." % region)
            return True

        if not super().add_parent(region, index=index, after=after,
                                  before=before, set_changed=False):
            return False

        index_shift = None
        if update_indices:
            # Determine index, where insertion took place.
            index = 0
            for parent in self.parents:
                if parent is region:
                    break
                index += parent.datapoints
    
            # Create index_shift: (caller, name, index, shift)
            # Virtually, adding a region, is as if one would shift the stop of the
            # preceding region by region.datapoints.
            index_shift = (self, 'stop', index, region.datapoints)

        # Inform children about index_shift
        self.set_changed(level=1, index_shift=index_shift)

        return True

    def remove_parent(self, region, update_indices=True):
        if region is None:
            return False

        # Possible causes of an index shift:
        #    1. index shift of parent
        #    2. adding or removing a parent
        # This function is called in case 2. For handling of case 1 see method
        # self.member_changed().

        # Determine index, where removal will take place.
        index = 0
        for parent in self.parents:
            if parent is region:
                break
            index += parent.datapoints

        # Remove parent
        # super().remove_parent() would call set_changed() itself, without any
        # parameters. This  means, if a parent was removed, set_changed() would
        # trigger a recalculation of Modification, regardless of whether the
        # change affected the region selected by the child View, or not.
        # Therefore, disable generic set_changed() mechanism and initiate a
        # specific one with the index_shift information.
        if not super().remove_parent(region, set_changed=False):
            return False

        index_shift = None
        if update_indices:
            # Create index shift
            # Virtually, removing a region is as if one would shift the start of a
            # following region by region.datapoints.
            index_shift = (self, 'start', index, region.datapoints)

        # Inform children about index_shift,
        self.set_changed(level=1, index_shift=index_shift)

        return True

    def _get_data_uncached(self, samples, traces_idx, copy=True):
        """
        Return concatenated data from all parents.
        """

        # determine num_traces
        if isinstance(traces_idx, slice):
            stop = traces_idx.stop
            if stop < 0:
                stop = stop + len(self.traces)
            num_traces = int(np.ceil((stop - traces_idx.start)
                                     / traces_idx.step))
        else:
            # traces is an np.array
            num_traces = len(traces_idx)

        # determine the global boundaries of requested samples and datapoints
        if isinstance(samples, slice):
            s_start = samples.start
            s_stop = samples.stop
            s_step = samples.step
            datapoints = int(np.ceil((s_stop - s_start) / s_step))
        else:
            # samples is an np.array
            s_start = samples.min()
            s_stop = samples.max() + 1
            # datapoints = len(samples)
            # datapoints from every single parent need to be appended on the
            # fly (np.vstack)
            datapoints = 0

        # Go through all parents and get data from all parents, whose
        # datapoints fall within the requested samples region.
        # Initialize running indices start and stop, which keep track of
        # samples served by every single parent.
        start = 0
        stop = 0
        # initialize array, which will hold all data from parents
        data = np.empty((datapoints, num_traces))
        for parent in self.parents:
            # set lower bound of index of data served by current parent
            # -> start of current parent is stop of previous parent
            start = stop
            # set upper bound of index of data served by current parent
            # -> stop of current parent is:
            stop = start + parent.datapoints
            # Check, whether data of current parent is requested by samples
            if start < s_stop and stop > s_start:
                # create sub fraction of samples for current parent
                if isinstance(samples, slice):
                    # treat steps/decimate
                    if start <= s_start:
                        # requested data served by the first parent -> s_start
                        # is index of first datapoint -> no shift, correct
                        # sub_start by start
                        req_start = 0
                        sub_start = s_start - start
                    else:
                        # requested data served by a subsequent parent -> start
                        # is index of first datapoint -> shift start of data
                        # according to steps left over by num_samples % s_step
                        #  either some steps left over -> shift is ->
                        #   shift = s_step - steps_left_over
                        #  or if no stpes left over -> no shift ->
                        #   shift = s_step - s_step
                        shift = s_step - ((start - s_start) % s_step or s_step)
                        req_start = (start - s_start + shift) / s_step
                        sub_start = shift
                    req_stop = np.ceil((min(stop, s_stop) - s_start) / s_step)
                    sub_stop = min(stop, s_stop) - start
                    requested = slice(int(req_start), int(req_stop))
                    sub_samples = slice(sub_start, sub_stop, s_step)
                    # get data from parent
                    data[requested] = parent._get_data(sub_samples, traces_idx,
                                                       copy)
                else:
                    # samples is an np.array:
                    #   minimum index can be equal to start (0) or greater
                    #   maximum index has to be lesser than stop (datapoints)
                    requested = np.logical_and(samples >= start,
                                               samples < stop)
                    sub_samples = samples[requested] - start
                    # get data from parent and append to data
                    p_data = parent._get_data(sub_samples, traces_idx, copy)
                    data = np.vstack((data, p_data))

        return data
