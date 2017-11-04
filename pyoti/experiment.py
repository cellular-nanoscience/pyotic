# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:50:29 2014

@author: Tobias Jachowski
"""
import collections
import getpass
import os
import shutil
import tempfile
import transaction
from functools import wraps, update_wrapper
from zc.lockfile import LockError
from ZODB import FileStorage, DB
from ZODB.POSException import ConnectionStateError

from . import config as cf
from . import helpers as hp
from . import update
from . import version
from .calibration import calibration as cb
from .calibration import Calibration, CalibrationSource
from .graph import GraphMember
from .gui import GRS
from .modification import Modification
from .picklable import Attributes
from .region import record as rc
from .region import Region, View, Record
from pyoti.plugins.modifications.generic import GenericMod


def if_open(func):
    @wraps(func)
    def func_wrapper(self, *args, **kwargs):
        if self.is_open:
            return func(self, *args, **kwargs)
        print("Experiment is closed!")
    return func_wrapper


def if_closed_return(value_or_func=None):
    def if_closed_return_decorator(func):
        @wraps(func)
        def func_wrapper(self, *args, **kwargs):
            if self.is_open:
                return func(self, *args, **kwargs)
            print("Experiment is closed!")
            try:
                value = value_or_func()
            except:
                value = value_or_func
            return value
        return func_wrapper
    return if_closed_return_decorator


class Experiment(object):
    """
    Provide functions to analyse tweezers data and save the analysis in a file:

    * create records (i.e. reading tweezers data),
    * construct groups, consisting of views and/or modifications:
        * create views (i.e. select timespan of data),
        * create modifications (i.e. modify the data),
        * change the timespan of views graphically,
        * change the parameters of modifications graphically,
        * concatenate and remove views (i.e. concatenate the data),
    * assemble/stack different groups,
    * manage records, views, and modifications,
    * delete records, views, and modifications.
    * load/save the records, views, and modifications and the graph structure

    Attributes
    ----------
    opened_last : Experiment
        References the `Experiment`, which opened the last experiment file.
    created_last : Experiment
        References the last created `Experiment`.
    """

    opened_last = None
    created_last = None

    def __init__(self, filename=None, directory=None, **kwargs):
        """
        Parameters
        ----------
        filename : str, optional
            The name of the experiment database file that should be opened.
            If filename is omitted, no file is opened.
        directory : str, optional
            The directory where the experiment is stored (defaults to '.').

        Attributes
        ----------
        filename : str
        version : str
        records : self.record(name)
            Offers easy access to Records via dotted notation:
            self.records.name is the same as calling self.record('name')
        views : self.view(name)
            Offers Easy access to Views via dotted notation: self.views.name is
            the same as calling self.view('name').
        modifications : self.modification(name)
            Offers Easy access to Modifications via dotted notation:
            self.modifications.name is the same as calling
            self.modification('name')
        cached_region : list of Regions
            The Region that is cached. It is stored in status.
        last_added_region : Region
            The Region that was added last to the experiment. It is stored in
            status.
        last_adjusted_region : Region
            The region that was adjusted last. It is stored in status.
        """
        self._grs = GRS()

        class AttributeCaller(object):
            def __init__(self, attribute_function, callback_function):
                self.attribute_function = attribute_function
                self.callback_function = callback_function
                update_wrapper(self, callback_function)

            def __getattr__(self, name):
                return self.attribute_function(name=name)

            def __call__(self, name=None, group=None):
                return self.callback_function(name=name, group=group)

        # Convenient objects to access records, views, and modifications via
        # dot notation: self.record.alpha or self.view.used, ...
        self.region = AttributeCaller(self._region, self._region)
        self.record = AttributeCaller(self._record, self._record)
        self.view = AttributeCaller(self._view, self._view)
        self.modification = AttributeCaller(self._modification,
                                            self._modification)
        # Backwards compatibility for "plural" notation
        self.regions = AttributeCaller(self._region, self._regions)
        self.records = AttributeCaller(self._record, self._records)
        self.views = AttributeCaller(self._view, self._views)
        self.modifications = AttributeCaller(self._modification,
                                             self._modifications)

        # Initialise attributes
        self._filename = None

        self._storage = None
        self._db = None
        self._transaction_manager = transaction.TransactionManager()
        self._connection = None

        self._tempdir = None

        if filename is not None:
            self.open(filename, directory)

        type(self).created_last = self

    def open(self, filename=None, directory=None, verbose=True):
        """
        Open experiment database file.
        """
        if self.is_open:
            if verbose:
                print("The following experiment file is already opened:\n",
                      "  '%s'\n" % self.filename,
                      "Close it, before you open a new file.")
            return False

        # Get and set database filename
        filename, absdir, absfile = hp.file_and_dir(filename, directory)

        if absfile is None:
            # Create a tempfile
            self._tempdir = tempfile.mkdtemp(prefix='pyoti-',)
            self._filename = tempfile.mktemp(suffix='.fs', prefix='',
                                             dir=self._tempdir)
            if verbose:
                print("No filename given. Switch on temporary mode.")
        else:
            self._filename = absfile

        if verbose:
            # Print whether a new experiment file was created or an old one was
            # opened
            if os.path.isfile(self.filename):
                print("Open experiment file:\n",
                      "  '%s'" % self.filename)
            else:
                print("Create experiment file:\n",
                      "  '%s'" % self.filename)

        # Try to open/create the experiment file
        try:
            # create/open db and TransactionManager
            self._storage = FileStorage.FileStorage(self.filename)
            self._db = DB(self._storage)
            self._connection = self._db.open(
                transaction_manager=self._transaction_manager)
        except LockError as err:
            if verbose:
                print("Cannot open the experiment file!\n",
                      "The database is locked by another process or user.\n",
                      "The original error message is: '%s'" % (err))
            return False
        except FileNotFoundError as err:
            if verbose:
                print("Cannot open the experiment file!\n",
                      "The file cannot be created!\n",
                      "The original error message is: '%s'" % (err))
            return False
        except Exception:
            raise  # raise unknown errors

        # Initialise the database of a new experiment file
        self._init_db()

        # Update the temporary absdir of graphroot. Make the attribute
        # volatile, so that ZODB doesn't store it.
        filename, absdir, absfile = hp.file_and_dir(self.filename)
        self._graphroot._v_absdir = absdir

        # Update the database of an experiment file created with an old version
        update.update_db(self)

        # Save the current timestamp
        # now = datetime.datetime.now()
        # current_date = now.ctime()
        # current_date_iso = now.isoformat()
        # global _current
        type(self).opened_last = self
        return True

    @if_open
    def _init_db(self, autosave=True):
        # create status, graphroot and version if new experiment
        if ('version' not in self._dbroot
            and ('views' in self._dbroot
                 or 'modifications' in self._dbroot
                 or 'records' in self._dbroot)):
            # Saving the database before updating it, would persistently break
            # the tweezersdata objects.
            # self._dbroot['version'] = '0.0.1'
            # autosave = False
            self.close()
            raise TypeError(("The experiment file '%s' was created with the "
                             "tweezersdata package.\n"
                             "This version existed before the pyoti, "
                             "investigator, or tethered_bead package.\n"
                             "Please, load the experiment file with the "
                             "tweezersdata package.") % self.filename)
        if 'status' not in self._dbroot:
            self._dbroot['status'] = Attributes()
            self._status.cached_region = []
            self._status.last_added_region = None
            self._status.last_adjusted_view = None
            autosave and self._save(description="Add status")
        if 'graphroot' not in self._dbroot:
            self._dbroot['graphroot'] = GraphMember(max_parents=0,
                                                    name='graphroot',
                                                    group='graphroot')
            autosave and self._save(description="Add graphroot")
        if 'version' not in self._dbroot:
            # Keep track of sourcecode used to create the database
            self._dbroot['version'] = version()
            autosave and self._save(description="Add version")

    @if_closed_return(True)
    def close(self, verbose=True, discard_temp=True):
        """
        Close the file the Experiment is stored in.
        """
        # Check if database was opened in temporary mode. If so, abort, close,
        # and delete temporary files.
        if discard_temp and self._tempdir is not None:
            print("Temporary experiment file to be deleted:\n  '%s'"
                  % self.filename)
            self.abort()
            self._close(verbose=verbose)
            self._rmtemp()
            return True

        return self._close(verbose=verbose)

    @if_closed_return(True)
    def _close(self, verbose=True):
        try:
            self._connection.close()
            self._db.close()
            self._storage.close()
            self._connection = None
            self._transaction_manager = transaction.TransactionManager()
            self._db = None
            self._storage = None
            self._filename = None
        except ConnectionStateError:
            if verbose:
                print("Cannot close experiment file due to unsaved changes.\n",
                      "If you really want to close without saving, execute",
                      "abort() first.")
            return False
        except Exception:
            raise  # raise unknown errors
        else:
            if verbose:
                print("Experiment closed.")
        return True

    def _rmtemp(self):
        if self._tempdir is None:
            return
        print("Switching off temporary mode.")
        shutil.rmtree(self._tempdir, ignore_errors=True)
        self._tempdir = None

    @if_open
    def create_record(self, name='alpha', group=None, calibration=None,
                      cfgfile=None, auto_add=True, **kwargs):
        """
        Open existing or create a new `Record`.

        Parameters
        ----------
        name : str
            The name of the `Record`. If a `Record` with the given `name`
            already exists in this `Experiment`, no new `Record` is created and
            the existing `Record` is returned.
        group : str, optional
            The name of the group, Defaults to `name`.
        calibration : Calibration or CalibrationSource, optional
            The `calibration` is used to create the new `Record`. If you
            provide no calibration, a generic calibration will be automatically
            created.
        cfgfile : str, optional
            The configuration file describing the setup specific parameters.
            Also see method `region.record.create_record()`.
        auto_add : bool, optional
            If `auto_add` == True, the new `Record` is automatically added to
            this `Experiment` by calling `add_record()`.
        **kwargs
            Keyword arguments are passed to function
            `pyoti.region.record.create_record()`

        Returns
        -------
        Record

        Notes
        -----
        See attribute `records` and method `record()` for retrieving
        already added records.

        """
        # set default group to name
        group = group or name

        # Try to find an already created record
        record = self.record(name=name, group=group)

        # create a new record with name and group
        if record is None:
            # Create a standard calibration
            if calibration is None:
                calibration = create_calibration()

            if isinstance(calibration, CalibrationSource):
                calibration = Calibration(calibration)

            # create the record
            record = rc.create_record(calibration, name=name, group=group,
                                      cfgfile=cfgfile, root=self._graphroot,
                                      **kwargs)

        # record could not be created
        if record is None:
            print("Could not create record")
            return None

        # automatically add the record to this experiment (`self._graphroot`)
        if auto_add:
            record = self.add_record(record=record)

        return record

    @if_open
    def add_record(self, record=None, set_last_added_region=True):
        """
        Add a Record to self, if a Record with the same name does not exist.

        Parameters
        ----------
        record : Record
            The Record to be added to self. Do not add the `record`, if a
            Record with the same name is already saved in self. Instead, return
            the already saved Record.
        set_last_added_region : bool
            If True and `record` was added, set self.last_added_region to the
            added Record.

        Returns
        -------
        Record
            This is either the added or already exising Record. If `record`
            could not be added nor a Record with the same name already exists,
            return None.
        """
        # read in record
        # If record already exists, do nothing.

        if record is None or not isinstance(record, Record):
            print("Can only add instances of `Record`!")
            return None

        # Check if record with same name already exists
        record_ = self.record(name=record.name, group=record.group)

        if record_ is None:
            # record with that name does not exist yet -> add it to self
            print("Adding record '" + record.name + "':")
            # add record to root graph member
            self._graphroot.add_child(record, set_changed=False,
                                      bidirectional=False)
            if set_last_added_region:
                self.last_added_region = record

            # Set datasource._root to self._graphroot for foreign records
            if (hasattr(record.datasource, '_root')
               and record.datasource._root is None):
                record.datasource._root = self._graphroot
        else:
            # record with that name does exist already -> set record to the
            # existing one
            record = record_
            print("Using already existing record '%s' and existing "
                  "calibration:" % (record.name))

        print(("  '" + record.datasource.name + "'"))
        print(("  '" + record.calibration.calibsource.name + "'"))

        return record

    @if_open
    def append_group(self, name, group=None, group_type='selection',
                     adjust=True, cfgfile=None, **kwargs):
        """
        Append a group to the last added group of this experiment. For details,
        parameters, and returns see method `add_group()`.
        """
        parent_region = self.last_added_region
        self.add_group(name,
                       parent_region,
                       group=group,
                       group_type=group_type,
                       adjust=adjust,
                       cfgfile=cfgfile,
                       **kwargs)

    @if_open
    def add_group(self,
                  name,
                  parent_region,
                  group=None,
                  parent_group=None,
                  group_type='selection',
                  adjust=True,
                  cfgfile=None,
                  traces=None,
                  description=None,
                  modclass=None,
                  **kwargs):
        """
        Add a group which will retrieve its data from a given `parent`.
        The type of the group is determined by the parameter `group_type`.
        The different types of groups are defined in the configfile (`cfgfile`,
        defaults to "groups.cfg").

        A group can be either (a) created for solely selecting a timespan of
        the data retrieved by the parent or (b) additionally modifying the
        data. In the case (a) `add_group()` creates a View by calling the
        method `add_view()`. In the case (b)  `add_group()` creates a View and
        additionally a Modification, by calling the method
        `add_modification()`.


        Parameters
        ----------
        name : str
            The name of the members of the group. For Modifications and their
            Views based, the name will be changed to name + "_mod".
        parent_region : str or GraphMember
            Name of the parent region, a GraphMember with the name, or the
            parent Region itself.
        group : str, optional
            The name of the group, Defaults to `name`.
        parent_group : str, optional
            Name of the group of the `parent_region`. Defaults to
            None or `parent_region.group`.
        group_type : str
            The name of the type of group defined in the file `cfgfile`.
        adjust : bool
            If True, the user is presented a graphical selection dialog (GRS)
            to select the timespan of the new group (i.e. the View).
            See also method `adjust_view()` and `adjust_modification()`.
        cfgfile : str, filepath, optional
            The configfile where the different types of groups are defined.
            Values can be overwritten by parameters `traces`, `description`,
            and `modclass`
            Defaults to "groups.cfg"
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the new group (i.e. View).
            Takes precedence over the settings of `cfgfile`. Defaults to all
            traces available in `parent_region`.
        description : str
            Short description of what the selection of the timespan is used
            for. Takes precedence over the settings of `cfgfile`.
        modclass : Modification
            If a class for a Modification is given, create a modifying group
            by calling the method `add_modifiation()`, else create a selecting
            group by calling the method `add_view()`. Defaults to None.
            Takes precedence over the settings of `cfgfile`.
        **kwargs : optional
            Additional Arguments that are passed either to the method
            `add_view()` or the method `add_modification()`.

        Returns
        -------
        View or Modification
            Depending on the type of group, return either the added `View` or
            the `Modification`.
        """
        cfgfile = cfgfile or 'groups.cfg'
        cfg = cf.read_cfg_file(cfgfile)

        description = description or cf.get_cfg_option(cfg, sec=group_type,
                                                       opt='description')
        traces = (traces or cf.get_cfg_list(cfg, sec=group_type, opt='traces')
                  or None)
        modclass = modclass or cf.get_cfg_class(cfg, sec=group_type,
                                    std_mod='.plugins.modifications.generic',
                                    cls_opt='modclass')

        if not modclass:
            return self.add_view(name, parent_region, group=group,
                                 parent_group=parent_group, adjust=adjust,
                                 traces=traces, description=description,
                                 **kwargs)
        else:
            return self.add_modification(name, group=group,
                                         modclass=modclass,
                                         parent_region=parent_region,
                                         parent_group=parent_group,
                                         adjust=adjust,
                                         traces=traces,
                                         description=description, **kwargs)

    @if_open
    def add_view(self, name, parent_region, group=None, parent_group=None,
                 adjust=True, traces=None, description=None, **kwargs):
        """
        Create and adjust a View with `parent_region` as a parent. The View
        will retrieve its data from the `parent_region`. The timespan is set to
        the same as the one from the parent. The new `View` is automatically
        added to this `Experiment`. If a `View` with the `name` already exists
        in this `Experiment`, do not create a new `View`, but instead retrieve
        and return the existing one. Let the user choose the timespan of the
        View, if adjust is True.

        Parameters
        ----------
        name : str
            The name of the View.
        parent_region : str or GraphMember
            Name of the parent region, a GraphMember with the name, or the
            parent Region itself.
        group : str, optional
            Name of the group, defaults to `name`.
        parent_group : str, optional
            Name of the group of the `parent_region`. Defaults to
            None or `parent_region.group`.
        adjust : bool
            If True, the user is presented a graphical selection dialog (GRS)
            to select the timespan of the new View.
            See also method `adjust_view()`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the new View. Defaults to all traces available in the
            Region determined by `parent_region`.
        description : str
            Short description of what the selection of the timespan is used
            for. Takes precedence over the settings of `cfgfile`.
        **kwargs : optional
            Additional arguments to initialize a View. See `View' for
            details. Only use, if you really now what you ar doing.

        Returns
        -------
        View
        """
        group = group or name
        view = self._add_view(name, group, parent_region,
                              parent_group=parent_group, **kwargs)
        if adjust:
            self.adjust_view(view, group=group, traces=traces,
                             description=description)

        return view

    @if_open
    def _add_view(self, name, group, parent_region, parent_group=None,
                  set_last_added_region=True, **kwargs):
        """
        Create a View with `parent_region` as a parent. The View will retrieve
        its data from the `parent_region`. The timespan is set to the same as
        the one from the parent. The new `View` is automatically added to this
        `Experiment`. If a `View` with the `name` already exists in this
        `Experiment`, do not create a new `View`, but instead retrieve and
        return the existing one.

        Parameters
        ----------
        name : str
            The name of the View.
        group : str
            The name of the group the View should be part of.
        parent_region : str or GraphMember or collections.Iterable of these
            Name of the parent region, a GraphMember with the name, or the
            parent Region itself.
        parent_group : str, optional
            Name of the group of the parent_region. Defaults to None or
            `parent_region.group`.
        set_last_added_region : bool
            If True and a new View was added, set self.last_added_region to the
            added View.
        **kwargs : optional
            Additional arguments to initialize a View. See `View' for
            details. Only use, if you really now what you ar doing.

        Returns
        -------
        View
        """
        # try to get already existing view
        view = self.view(name=name, group=group)

        # view does not exist -> create new view
        if view is None:
            parent = self.region(name=parent_region, group=parent_group)
            if parent is not None:
                # switch on caching of this Region's parent
                # and switch off caching of previous Region's parent
                self.set_cached_region(parent)

                print(("Creating view '" + name + "'."))
                view = View(parent=parent, caching=False, name=name,
                            group=group, **kwargs)
                if set_last_added_region:
                    self.last_added_region = view
            else:
                print("Could not find parent %s" % (parent_region))
        else:
            print(("Using already existing view '" + name + "':"))

        return view

    @if_open
    def show_view(self, view, group=None, traces=None, description=None,
                  xlim=None):
        """
        Graphically inspect a View.

        Parameters
        ----------
        view : str or GraphMember
            Name of the View, a GraphMember with the name, or the View itself.
        group : str, optional
            The name of the group. Defaults to None or `view.group`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the View. Defaults to all traces available in the
            View determined by `view`.
        description : str
            Short description of what the selection of the timespan is used
            for.
        xlim : (float, float), optional
            The timespan (s) plotted.
            Defaults to (View.tmin, View.tmax) of the View determined by
            `view`.
        """
        view = self.view(name=view, group=group)
        if view is not None:
            # if traces is None or (isinstance(traces, list)
            # and len(traces) == 0):
            traces = view.traces_available(traces)
            # tc.normalize(traces)
            # traces = traces or view.traces
            if description is None:
                description = view.name
            self._grs.init_ifig(view.timevector,
                                view.get_data(traces=traces, copy=False),
                                view.samplingrate,
                                traces,
                                xlim=xlim,
                                description=description,
                                select=False)

    @if_open
    def adjust_view(self, view, group=None, traces=None, description=None,
                    xlim=None, set_last_adjusted_view=True):
        """
        Graphically adjust the timespan of a View.

        Parameters
        ----------
        view : str or GraphMember
            Name of the View, a GraphMember with the name, or the View itself.
        group : str, optional
            The name of the group. Defaults to None or `view.group`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the View. Defaults to all traces available in the
            View determined by `view`.
        description : str
            Short description of what the selection of the timespan is used
            for.
        xlim : (float, float), optional
            The timespan (s) plotted the user can choose the timespan from.
            Defaults to (View.parent.tmin, View.parent.tmax) of the View
            determined by `view`.
        set_last_adjusted_view : bool
            If True, set self.last_adjusted_view to the View determined by
            `view`.
        """
        view = self.view(name=view, group=group)
        if view is not None:
            # if traces is None or (isinstance(traces, list)
            # and len(traces) == 0):
            traces = view.traces_available(traces)
            # tc.normalize(traces)
            # traces = traces or view.traces
            if description is None:
                description = view.name
            self._grs.init_ifig(view.parent.timevector,
                                view.parent.get_data(traces=traces,
                                                     copy=False),
                                view.parent.samplingrate,
                                traces,
                                tmin=view.tmin,
                                tmax=view.tmax,
                                xlim=xlim,
                                onselect_cb=view.set_timespan,
                                description=description)

            if set_last_adjusted_view:
                self.last_adjusted_view = view

    @if_open
    def add_modification(self, name, group=None, modclass=None,
                         parent_region=None, parent_group=None,
                         view_apply=None, mod_name_append='_mod', adjust=True,
                         traces=None, description=None, **kwargs):
        """
        Create a `Modification` and associate it to a View the Modification is
        applied to (View applied) and a View whose data is used to calculate
        its parameters from (View based).
        The name of the View applied is `name`, or `apply_view`, if given. The
        name of the View based and the Modification is determined by
        concatenating `name` and `mod_name_append`.
        If either of the Views does not exist, automatically create them with
        `parent_region` as the parent Region (details see `add_view()`)
        If `adjust`=True, let the user graphically adjust the timespan of the
        View based on. Otherwise, leave the timespan untouched, if the View
        based already exists, or set the timespan to the same as the
        `parent_region` if the View based has to be created.

        Parameters
        ----------
        name : str
            The first part of the name of the Modification, and, if
            `apply_view` is not given, the name of the `View` the Modification
            should be applied to. If a View with `name` name does not exist,
            create a new one. See `add_view()` for details.
        group : str, optional
            The name of the group. Defaults to `name`.
        modclass : Modification
            The class of the Modification that should be used.
        parent_region : str or GraphMember, optional
            Name of the parent region, a GraphMember with the name, or the
            parent Region itself. `parent_region` will be used for the creation
            of the View applied as well as the View based. If it is not
            supplied, the parent(s) of the View applied are used. See
            `add_view()` for details.
        parent_group : str, optional
            Name of the group of the parent_region. Defaults to None or
            `parent_region.group`.
        view_apply : str or GraphMember, optional
            The name of the `View` the Modification should be applied to. If a
            View with name `name` does not exist, create a new one, with
            parent `parent_region`. See `add_view()` for details.
        mod_name_append : str
            The name of the Modification and the View based are determined by
            concatenating 'name' and 'mod_name_append'.
        adjust : bool
            If True, the user is presented a graphical selection dialog (GRS)
            to select the timespan of the View based. See also method
            `adjust_modification()`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the View. Defaults to all traces available in the
            Region determined by `parent_region`.
        description : str
            Short description of what the selection of the timespan is used
            for.
        **kwargs
            Keyword arguments are passed on to initializing class `modclass`
        """
        # Create or get a view the modification is applied to
        if view_apply is None and parent_region is None:
            raise TypeError("add_modification() is missing keyword argument "
                            "`parent_region` or `view_apply`")

        group = group or name

        # Create a new view_apply, parent_region needs to be given
        if view_apply is None:
            view_apply = self._add_view(name, group, parent_region,
                                        parent_group=parent_group)

        # Set parent region(s) for the view_based of modification to be the
        # same as the ones from view_apply, view_apply needs to be given
        if parent_region is None:
            # Get existing view
            view_apply = self.view(name=view_apply, group=group)
            parent_region = view_apply.parent.parents

        # Create modification
        mod_name = name + mod_name_append

        # Set modclass to be defaulting to GenericMod
        modclass = modclass or GenericMod

        modification = self._add_modification(mod_name, group, modclass,
                                              view_apply, parent_region,
                                              parent_group=parent_group,
                                              **kwargs)
        if adjust:
            self.adjust_modification(modification, group=group, traces=traces,
                                     description=description)

        return modification

    @if_open
    def _add_modification(self, name, group, modclass, view_apply,
                          parent_region, parent_group=None, **kwargs):
        """
        Create a `Modification` and associate it to a View the Modification is
        applied to (View applied) and a View whose data is used to calculate
        its parameters from (View based).
        The name of the Modification and, if `parent_region` is given, the View
        based is `name`. The class used to create the Modification is
        `modclass`.
        If the Modification with the name `name` does exist, simply return the
        modification.

        Parameters
        ----------
        name : str
            The name of the Modification and the View based. If a Modification
            with name `name` does not exist, create one. If `parent_region` is
            given, use the name also for the View based, that is created, too,
            if a View with the same name does not exist yet.
        group : str
            Name of the group.
        modclass : Modification
            The class of the Modification that should be used.
        view_apply : str or GraphMember
            The name of the `View` applied.
        parent_region : str or GraphMember or Iterable of str/GraphMember
            The name of the parent Region(s) that is used as a parent for View
            based, if Modification and View based do not exist and have to be
            created.
        parent_group : str, optional
            Name of the group of the parent_region. Defaults to None or
            `parent_region.group`.
        **kwargs
            Keyword arguments are passed on to initializing of class
            `modclass`, when a new Modification has to be created.

        Returns
        -------
        Modification
            The created or found Modification.
        """
        modification = self.modification(name=name, group=group)

        if modification is None:
            # Create new modificiation
            view_apply = self.view(name=view_apply, group=group)
            parent_region = self.region(name=parent_region, group=parent_group)

            if view_apply is not None and parent_region is not None:
                # modification should  have a parent_region/view_based so that
                # it can be referenced by traversing the graph from top to
                # bottom

                # switch on caching of this Region's parent and switch off
                # caching of previous Region's parent
                self.set_cached_region(parent_region)

                view_based = self._add_view(name, group, parent_region,
                                            set_last_added_region=False)

                print(("Creating modification '" + name + "'."))
                # Instantiate a Modification subclass
                # e.g. 'Offset' would correspond to subclass Offset
                modification = modclass(view_based=view_based,
                                        view_apply=view_apply, name=name,
                                        group=group, **kwargs)
            else:
                print("Could not find/use view_apply %s and/or parent_region "
                      "%s." % (view_apply, parent_region))
        else:
            print(("Using already existing modification '" + name + "':"))

        # initialize the modification of the modification
        # if this would not be done directly after having added th
        # modification the modification would first be calculated, when trying
        # to get data from the View the modification is applied to.
        # modification.evaluate()

        return modification

    @if_open
    def adjust_modification(self, modification, group=None, traces=None,
                            description=None, xlim=None):
        """
        Adjust the timespan of a View based of a Modification.

        Parameters
        ----------
        modification : str or GraphMember
            Name of the Modfication whose View based should be adjusted.
        group : str, optional
            The name of the group. Defaults to None or `modification.group`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the View. Defaults to all traces available in the
            View determined by `view`.
        description : str
            Short description of what the selection of the timespan is used
            for.
        xlim : (float, float), optional
            The timespan (s) plotted the user can choose the timespan from.
            Defaults to (View.parent.tmin, View.parent.tmax) of the View
            determined by `view`.
        """
        modification = self.modification(modification, group=group)

        self.adjust_view(modification.view_based, group=group, traces=traces,
                         description=description, xlim=xlim)

        try:
            # Interact with the modification, after view has been selected
            modification.interact()
        except Exception as exc:
            print("Something went wrong!\n",
                  "Can not interact with the modification!\n",
                  "The original error message is: '%s'" % (exc))

    @if_open
    def concatenate(self, name, first_region, second_region, group=None,
                    first_group=None, second_group=None, adjust=True,
                    traces=None, description=None):
        """
        Concatenate two Regions as parent Regions for a new View with name
        `name`. If the View with the name `name` already exists, append the
        `second_region` to the parent Regions of the View with the name `name`.

        Parameters
        ----------
        name : str
            Name of the View, which should have the data of the concatenated
            Regions.
        first_region : str or GraphMember
            Name of the Region, a GraphMember with the name, or the Region
            itself, whose data should be concatenated with second_region.
        second_region : str or GraphMember
            Name of the Region, a GraphMember with the name, or the Region
            itself, whose data should be concatenated with first_region.
        group : str, optional
            Name of the group for the new concatenated view. Defaults to
            `name`.
        first_group : str, optional
            Name of the group of the first_region. Defaults to None or
            `first_region.group`.
        second_group : str, optional
            Name of the group of the second_region. Defaults to None or
            `second_region.group`.
        adjust : bool
            If True, the user is presented a graphical selection dialog (GRS)
            to select the timespan of the new View.
            See also method `adjust_view()`.
        traces : str, int, list of str/int, or slice
            Select the traces the user is shown to adjust the timespan from.
            The selection of the traces does *not* influence the traces
            contained in the View. Defaults to all traces available in the
            View determined by `view`.
        description : str
            Short description of what the selection of the timespan is used
            for.

        Returns
        -------
        View
            The View which presents the concatenated data.
        """
        group = group or name
        view = self._add_view(name, group, first_region,
                              parent_group=first_group)
        self.append_to(view, second_region, group=group,
                       append_group=second_group)
        if adjust:
            self.adjust_view(view, group=group, traces=traces,
                             description=description)
        return view

    @if_open
    def append_to(self, name, append_region, group=None, append_group=None,
                  after=None, after_group=None, before=None,
                  before_group=None, update_indices=True):
        """
        Appends a Region as a parent to a View.

        Parameters
        ----------
        name : str
            Name of a View a new parent should be appended to.
        append_region : str or GraphMember
            Name of a Region, a GraphMember with the name, or the Region
            itself, which should be appended as a parent.
        group : str, optional
            Name of the group for the view the data should be appended to.
            Defaults to `name`.
        append_group : str, optional
            Name of the group of the append_region. Defaults to None or
            `append_region.group`.
        after : str or GraphMember
            Name of a Region, a GraphMember with the name, or the Region
            itself, after which the new Region should be appended to.
        after_group : str, optional
            Name of the group of the after region. Defaults to None or
            `after_region.group`.
        before : str or GraphMember
            Name of a Region, a GraphMember with the name, or the Region
            itself, before which the new Region should be appended to. If after
            was given and could be found, before is ignored.
        before_group : str, optional
            Name of the group of the before region. Defaults to None or
            `before_region.group`.
        update_indices : bool, optional
            Update the indices of the child Regions. Defaults to True.

        Returns
        -------
        bool
            True if the region was appended or is already a parent. Otherwise
            return False.
        """
        if group is None:
            if isinstance(name, GraphMember):
                group = name.group
            else:
                group = name
        view = self.view(name=name, group=group)
        append_region = self.region(name=append_region, group=append_group)
        if view is not None and append_region is not None:
            after = self.region(name=after, group=after_group)
            before = self.region(name=before, group=before_group)
            return view.parent.add_parent(append_region, after=after,
                                          before=before,
                                          update_indices=update_indices)
        return False

    @if_open
    def remove_record(self, record, group=None, dereference_children=True):
        """
        Remove a Record with the name `name`.

        Remove all references from record to its parent (self._graphroot) and
        vice versa.

        Be careful, if you remove a record, all children, which have the record
        as their only root, will be removed, too!

        Parameters
        ----------
        record : str or GraphMember
            Name of the Record, a GraphMember with the name, or the Record
            itself that should be removed.
        group : str, optional
            Name of the group of the record. Defaults to None or
            `record.group`.

        Returns
        -------
        Record
            The removed record.
        """
        # removing a record
        record = self.record(name=record, group=group)
        if record is not None:
            # remove record from root graph member
            self._graphroot.remove_child(record)
            # dereference children from record
            if dereference_children:
                for child in record.children:
                    record.remove_child(child)
            # Remove reference of datasource to self._graphroot
            if hasattr(record.datasource, '_root'):
                record.datasource._root = None
        return record

    @if_open
    def remove_view(self, view, group=None):
        """
        Remove a View or a View based and its associated Modification with the
        name `view`.

        Remove all references from view to its parents (view.parent.parents
        and view.modifications_apply) and its children (view.children and
        view.modifications_based) and vice versa.

        Parameters
        ----------
        view : str or GraphMember
            Name of the View, a GraphMember with the name, or the View itself
            that should be removed.
        group : str, optional
            Name of the group of the view. Defaults to None or `view.group`.

        Returns
        -------
        View
            The removed view.

        Notes
        -----
        The structure for a View could look like this:

        parent.parents
              |
             view             <-- region
              |
          children
        """
        view = self.view(name=view, group=group)
        if view is not None:
            # First, delete relations to parents (view.parent is a MultiRegion)
            for parent in view.parent.parents:
                parent.remove_child(view.parent)
            for mod_apply in view.modifications_apply:
                mod_apply.remove_child(view)
            # Second, delete relations to children
            for mod_based in view.modifications_based:
                if mod_based.view_apply is not None:
                    mod_based.view_apply.remove_parent(mod_based)
            for child in view.children:
                child.remove_parent(view)

        return view

    @if_open
    def remove_modification(self, modification, group=None):
        """
        Remove a Modification with the name `modification`.

        Remove all references from modification to its parent (View based)
        and its child (View apply) and vice versa.

        Parameters
        ----------
        modification : str or GraphMember
            Name of the Modification, a GraphMember with the name, or the
            Modification itself that should be removed.
        group : str, optional
            Name of the group of the modification. Defaults to None or
            `modification.group`.

        Returns
        -------
        Modification
            The removed Modification.

        Notes
        -----
        The structure for a Modification could look like this:

        parent.parents
              |
          view_based          <-- region
              |  <- mod_based
              |
         modification
              |
              |  <- mod_apply
          view_apply
              |
          children
        """
        modification = self.modification(name=modification, group=group)
        if modification is not None:
            # First, delete relation to parent (view.parent is a MultiRegion)
            modification.view_based.remove_child(modification)
            # Second, delete relation to child
            modification.view_apply.remove_parent(modification)

        return modification

    @if_open
    def remove(self, name, group=None):
        """
        Call and return the return value of the function `remove_view()`.
        If a view with the name could not be found, try to remove a record.

        Parameters
        ----------
        name : str or GraphMember
        group : str, optional
            Name of the group of the graphmember to be removed. Defaults to
            None or `name.group`.
        """
        return self.remove_view(view=name, group=group) \
            or self.remove_record(record=name, group=group)

    @if_open
    def remove_from(self, view, remove_region, group=None, remove_group=None,
                    update_indices=True):
        """
        Remove a parent Region of a (concatenated) view. If the region is the
        last in the view, the view and all its descendants are removed, too.

        Removes all references from view to its parent remove_region and vice
        versa.

        Parameters
        ----------
        view : str or GraphMember
            Name of the View, a GraphMember with the name, or the View itself
            that the parent Region should be removed from.
        remove_region : str or GraphMember
            Name of the Region, a GraphMember with the name, or the Region
            itself that should be removed.
        group : str, optional
            Name of the group of the view. Defaults to None or `view.group`.
        remove_group : str, optional
            Name of the group of the remove region. Defaults to None or
            `remove_region.group`.
        update_indices : bool, optional
            Update the indices of the child Regions. Defaults to True.

        Returns
        -------
        Region
            The removed region.
        """
        view = self.view(name=view, group=group)
        remove_region = self.region(name=remove_region, group=remove_group)

        if view is not None and remove_region is not None:
            view.parent.remove_parent(remove_region,
                                      update_indices=update_indices)

        return remove_region

    @if_open
    def replace_in(self, view, remove_region, with_region, group=None,
                   remove_group=None, with_group=None, update_indices=True):
        """
        Replace a parent Region of a view with another one.

        Parameters
        ----------
        view : str or GraphMember
            Name of the View, a GraphMember with the name, or the View itself
            that the parent Region should be replaced in.
        remove_region
            Name of the Region, a GraphMember with the name, or the Region
            itself that should be replaced.
        with_region
            Name of the Region, a GraphMember with the name, or the Region
            itself that should replace the `remove_region`.
        group : str, optional
            Name of the group of the view. Defaults to None or `view.group`.
        remove_group : str, optional
            Name of the group of the remove region. Defaults to None or
            `remove_region.group`.
        with_group : str, optional
            Name of the group of the with region. Defaults to None or
            `with_region.group`.
        update_indices : bool, optional
            Update the indices of the child Regions. Defaults to True.

        Returns
        -------
        Region
            The replaced region.
        """
        view = self.view(name=view, group=group)
        with_region = self.region(name=with_region, group=with_group)
        remove_region = self.region(name=remove_region, group=remove_group)
        if (view is not None
                and remove_region is not None
                and with_region is not None
                and remove_region is not with_region):
            self.append_to(view, with_region, after=remove_region,
                           update_indices=update_indices)
            self.remove_from(view, remove_region,
                             update_indices=update_indices)

        return remove_region

    @if_open
    def members(self, name=None, group=None, instance_class=None, dft=False,
                level=-1):
        """
        Yields the members of group `group` with name `name` and class
        `instance_class`.

        Parameters
        ----------
        name : `instance_class`, GraphMember, or str, optional
            The name of the group that should be searched for in this
            Experiment.
        group : str, optional
            The name of the group, Defaults to None or `name.group`.
        instance_class : class, optional
            The class of the instance that is searched for and returned.
        dft : bool, optional
        level : int, optional

        Yields
        ------
        GraphMember
        """
        if name is not None \
                and (group is None
                     or (isinstance(name, GraphMember)
                         and name.group is group)) \
                and (instance_class is None
                     or isinstance(name, instance_class)):
            yield name
        else:
            if isinstance(name, GraphMember):
                name = name.name
            # get the members
            members = self._graphroot.members(name=name, group=group,
                                              instance_class=instance_class,
                                              includeself=False, dft=dft,
                                              level=level)
            for member in members:
                yield member

    @if_open
    def member(self, name=None, group=None, instance_class=None, dft=False,
               level=-1):
        """
        Returns the last member of group `group` with name `name` and class
        `instance_class`.

        Parameters
        ----------
        name : instance_class, GraphMember, or str, optional
            The name of the group that should be searched for in this
            Experiment.
        group : str, optional
            The name of the group, Defaults to None or `name.group`.
        instance_class : class, optional
            The class of the instance that is searched for and returned.
        dft : bool, optional
        level : int, optional

        Returns
        -------
        member : GraphMember
        """
        # First get the members
        members = self.members(name=name, group=group,
                               instance_class=instance_class, dft=dft,
                               level=level)
        # Usually, a group consists of the following parent -> child
        # combinations:
        # 1. Record
        # 2. MultiRegion -> View
        # 3. Modification -> View
        # 4. MultiRegion -> View -> Modfication -> View
        #
        # To always get the View and not the MultiRegion, if a Region is asked
        # for, return the last member found member
        last = None
        for member in members:
            last = member
        return last

    @if_open
    def _regions(self, name=None, group=None):
        """
        Finds and yields all Regions of this `Experiment`.

        Yields
        ------
        Region
        """
        return self.members(name=name, group=group, instance_class=Region)

    @if_open
    def _records(self, name=None, group=None):
        """
        Finds and yields all Records of this `Experiment`.

        Yields
        ------
        Record
        """
        return self.members(name=name, group=group, instance_class=Record,
                            level=1)

    @if_open
    def _views(self, name=None, group=None):
        """
        Finds and yields all Views of this `Experiment`.

        Yields
        ------
        View
        """
        return self.members(name=name, group=group, instance_class=View)

    @if_open
    def _modifications(self, name=None, group=None):
        """
        Finds and yields all Modifications of this `Experiment`.

        Yields
        ------
        Region
        """
        return self.members(name=name, group=group,
                            instance_class=Modification)

    @if_open
    def _region(self, name=None, group=None):
        """
        If a Region with the same name as `name` is already existing in self,
        return the already existing one. Otherwise, if `name` is a Region
        itself, return `name`. If None of this is True, return None.

        Parameters
        ----------
        name : str or GraphMember

        Returns
        -------
        Region or None

        Notes
        -----
        The returned Region can be either a Record or a View.
        """
        if name is None and group is None:
            return None
        return self.member(name=name, group=group, instance_class=Region)

    @if_open
    def _record(self, name=None, group=None):
        """
        If a Record with the same name as `name` is already existing in self,
        return the already existing one. Otherwise, if `name` is a Record
        itself, return `name`. If None of this is True, return None.

        Parameters
        ----------
        name : str or GraphMember

        Returns
        -------
        Record or None
        """
        if name is None and group is None:
            return None
        return self.member(name=name, group=group, instance_class=Record,
                           level=1)

    @if_open
    def _view(self, name=None, group=None):
        """
        If a View with the same name as `name` is already existing in self,
        return the already existing one. Otherwise, if `name` is a View
        itself, return `name`. If None of this is True, return None.

        Parameters
        ----------
        name : str or GraphMember

        Returns
        -------
        View or None
        """
        if name is None and group is None:
            return None
        return self.member(name=name, group=group, instance_class=View)

    @if_open
    def _modification(self, name=None, group=None):
        """
        If a Modification with the same name as `name` is already existing in
        self, return the already existing one. Otherwise, if `name` is a
        Modification itself, return `name`. If None of this is True, return
        None.

        Parameters
        ----------
        name : str or grap.GraphMember

        Returns
        -------
        Modification or None
        """
        if name is None and group is None:
            return None
        return self.member(name=name, group=group, instance_class=Modification)

    @if_open
    def print_status(self, notes=True, multiregion=False):
        """
        Prints a status of all Records, Views, and Modifications of this
        Experiment.
        """
        if notes and self.notes != "":
            print("Notes:")
            print(self.notes)
            print("-" * 119)
        print("name \t\t      group \t\tclass\t\t  caching/updated")
        print("-" * 119)
        if multiregion:
            for region in self.regions():
                print(str(region) + "\t"
                      + "c: " + str(region.caching))
        else:
            for record in self.records():
                print(str(record) + "\t"
                      + "c: " + str(record.caching))
            for view in self.views():
                print(str(view) + "\t"
                      + "c: " + str(view.caching))
        for modification in self.modifications():
            print(str(modification) + "\t"
                  + "u: " + str(modification.updated))

    @if_open
    def set_cached_region(self, region, group=None, update_cache=True):
        """
        Set the Region whose data should be cached and unset the caching of the
        previously cached Region (see attribute `self.cached_region`).

        Parameters
        ----------
        region : str or GraphMember or collections.Iterable of these
            Name of the Region, a GraphMember with the name, or the Region
            itself, whose caching should be switched on.
        group : str, optional
            The name of the group, Defaults to None or `region.group`.
        update_cache : bool
            Immidiately trigger an update of the cache of the new Region.
        """
        region = self.region(region, group=group)

        if not isinstance(region, collections.Iterable):
            region = [region]

        # switch on caching of this Region
        for rgn in region:
            rgn.caching = True
            if update_cache:
                rgn.update_cache()

        # switch off caching of previous Region
        if self.cached_region != region:
            for crgn in self.cached_region:
                if not isinstance(crgn, Record):
                    crgn.caching = False

        self.cached_region = region

    @if_open
    def has_group(self, group):
        """
        Check, if a group with the name `group` exist.

        A group can be either a View or a Modification with its View based.

        Parameters
        ----------
        group : str
            The name of the group.
        """
        return self._graphroot.group_root(group) is not None

    @if_open
    def save(self, pack=False, description=None, filename=None, directory=None,
             verbose=True):
        """
        Save the current status of the Experiment to file.

        Parameters
        ----------
        pack : bool, optional
            If pack is True, call `self.cleanup()`.
        description : str, optional
            Associate and stora a description along with all the changes from
            last till this saving point.
        filename : str, optional
            New filename to rename the experiment file to. Defaults to current
            filename. Only used if either `filename` or `directory` are set.
        directory : str, optional
            New directory to move the experiment file to. If tempfile, defaults
            to the current working directory. Otherwise, to the current
            directory of self.filename.
        """
        # Save the status of the experiment to file
        self._save(pack=pack, description=description)

        # Move the experiment file to a new location, if requested
        if filename is not None or directory is not None:
            self.move_to(filename=filename, directory=directory)

        # Warn user about the consequences, if temporary mode is switched on
        if self._tempdir is not None and verbose:
            print("ATTENTION: Experiment is in temporary mode, only. If you "
                  "close the experiment file, all work is lost! To switch off "
                  "the temporary mode and store the experiment file "
                  "permanently, either provide a `filename` and/or a "
                  "`directory`, or directly call the method `self.move_to` "
                  "prior to closing the experiment file.")

    @if_open
    def _save(self, pack=False, description=None):
        # Get current (or new) transaction
        txn = self._transaction_manager.get()

        # Set a description for the current transaction
        description = description or ""
        txn.note(description)
        # if sys.platform in ['linux', 'linux2', 'darwin']:
        #   doesn't work in notebook via internet
        #   user = os.getlogin()
        # else:
        # Set a user for the current transaction
        user = getpass.getuser()
        txn.setUser(user)

        # commit the transaction (i.e. save the experiment file)
        txn.commit()

        # Delete old undo data
        if pack:
            self.cleanup()

    @if_open
    def move_to(self, filename=None, directory=None):
        """
        Rename experiment file to a new name and/or move it to a now directory.

        Parameters
        ----------
        filename : str
            New filename to rename the experiment file to. Defaults to current
            filename.
        directory : str
            New directory to move the experiment file to. If tempfile, defaults
            to the current working directory. Otherwise, to the current
            directory of self.filename.

        Returns
        -------
        bool
        """
        # Get new and old database path
        old_filename, old_absdir, old_absfile = hp.file_and_dir(self.filename)
        if self._tempdir is None:
            # In non temp mode default directory to current one
            directory = directory or old_absdir
        filename, absdir, absfile = hp.file_and_dir(filename, directory)
        if filename is None:
            # If no filename was give, default to current filename
            filename, absdir, absfile = hp.file_and_dir(old_filename, absdir)

        if filename == old_filename and absdir == old_absdir:
            # No new filename and/or path, nothing to do, exit silently
            return False

        if not os.path.isdir(absdir):
            # directory where file should be stored does not exist
            print("A directory with the name '%s' does not exist." % absdir)
            return False

        if os.path.isfile(absfile):
            # another file with the same name already exists
            print("A file with the name '%s' already exists." % absfile)
            return False

        if os.path.isdir(absfile):
            # another directory with the same name already exists
            print("A directory with the name '%s' already exists." % absfile)
            return False

        # save current state of the experiment
        self._save()

        # Change ds._directory of all datasources to be relative to the new
        # absdir.
        changed = False
        for child in self._graphroot.children:
            if hasattr(child, 'datasource'):
                ds = child.datasource
                if ds._needs_file:
                    changed = True
                    ds._directory = os.path.relpath(ds.absdir, absdir)

        if changed:
            # save experiment with changed directories
            self._save(description="Change location of experiment to '%s'"
                       % absfile)

        # Close experiment file
        self.close(verbose=False, discard_temp=False)

        try:
            print("Move file from\n  '%s' to\n  '%s'"
                  % (old_absfile, absfile))
            # Move/rename experiment file to new location
            shutil.move(old_absfile, absfile)
        except:
            # File could not be moved, revert everything
            self.open(filename=old_absfile, verbose=False)
            if changed:
                self.undo(verbose=False)
            return False
        else:
            # Delete temporary ZODB files
            for ext in ['.lock', '.tmp', '.index', '.old']:
                try:
                    os.remove(''.join([old_absfile, ext]))
                except:
                    pass
            # delete tempdir if file was moved from tempdir
            self._rmtemp()

            # reopen the moved file
            self.open(filename=absfile, verbose=False)

        return True

    @if_open
    def abort(self):
        """
        Abort all changes of the Experiment back to the last saving point.
        """
        self._transaction_manager.abort()

    @if_open
    def cleanup(self):
        """
        Pack/clean the databasefile, the Experiment is stored in, to purge
        obsolete (undo) data.
        Changes of the database do not overwrite an existing object, they only
        create a new revision of that object. The benefit is, you may go back
        to an earlier revision, the drawback is, you have to pack the database
        from time to time to purge obsolete data.
        """
        self._db.pack()

    @if_open
    def undo(self, verbose=True):
        """
        Undo all changes of the Experiment back to the last saving point.
        """
        if not self._db.supportsUndo():
            if verbose:
                print("UNDO is not supported!")
            return

        # Get the last saved change
        undo_log = self._db.undoLog(0, -1)
        if len(undo_log) > 0:
            last_change = undo_log[0]
            if verbose:
                print("UNDOing: %s" % last_change)
            txn = self._transaction_manager.get()
            self._db.undo(last_change['id'], txn=txn)
            self._save(description="UNDO last change")
        else:
            if verbose:
                print('Nothing to UNDO.')

    @property
    def filename(self):
        """
        str
            The absolute path and filename of this experiment. It is either
            determined automatically upon initilization of the `Experiment`
            (temporary mode) or can alternatively be set by the parameters
            `filename` and `directory`.
        """
        return self._filename

    @property
    def is_open(self):
        """
        Check if experiment file is open.

        Returns
        -------
        bool
        """
        if '_connection' not in self.__dict__:
            return False
        return self._connection is not None

    @property
    @if_closed_return(True)
    def is_closed(self):
        """
        Check if experiment file is closed.

        Returns
        -------
        bool
        """
        return False

    @property
    @if_closed_return(dict)
    def _dbroot(self):
        return self._connection.root()

    @property
    @if_closed_return(version())
    def version(self):
        """
        str
            The version of the program (API and attributes) this experiment
            was created and saved.
        """
        return self._dbroot['version']

    @property
    @if_closed_return(dict)
    def _status(self):
        """
        picklable.Attributes
            self._status automatically stores any not yet existing attribute of
            the experiment. All stored attributes can be accessed via dotted
            notation and dict notation via self._status and dotted notation via
            experiment.

        Notes
        -----
        Make sure the attribute can be pickled and unpickled by ZODB/Python.
        """
        return self._dbroot['status']

    @property
    @if_closed_return()
    def notes(self):
        """
        Notes for this experiment.

        Returns
        -------
        str

        Examples
        --------
        To set the notes, assign the property `self.notes` as following:

        >>> self.notes = \"""
        >>> Put the multiline notes between the triple quotes.
        >>>
        >>> Empty lines are allowed.
        >>>
        >>> Empty lines at the beginning and the end of the text will be
        >>> stripped.
        >>> \"""
        """
        if 'notes' in self._status:
            return self._status.notes.strip()
        else:
            return ""

    @notes.setter
    @if_open
    def notes(self, notes):
        """
        Store notes for this experiment.

        Parameters
        ----------
        notes : str
            The notes to be stored
        """
        self._status.notes = notes

    @property
    @if_closed_return()
    def _graphroot(self):
        """
        GraphMember
            The root GraphMember, which holds the references (via a `Node`) to
            all Records of this experiment and subsequentially to all other
            Views and Modifications.
        """
        return self._dbroot['graphroot']

    def __getattr__(self, name):
        """
        Get the attributes that are stored in the `Attributes` status object
        of the Experiment.

        Parameters
        ----------
        name : str
            name is used as a key for self._status.
        """
        if self.is_open and name in self._status:
            return self._status[name]
        else:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Set an attribute that is stored in the `Attributes` status object of
        the Experiment.

        Parameters
        ----------
        name : str
            name is used as a key for self._status.
        value
            Object to be stored in self._status. Any picklable object will do.
        """
        # A __setattr__() hook will also override the Persistent __setattr__()
        # hook.  User code must treat it much like __getattribute__().  The
        # user-defined code must call _p_setattr() first to handle special
        # attributes; _p_setattr() takes the attribute name and value.  If it
        # returns True, Persistent handled the attribute. If not, the user code
        # can run. If the user code modifies the object’s state, it must be
        # assigned to _p_changed
        if (self.is_open  # Only set self._status.name variable if db is open
                          # 'self._connection' determines the open state.
                          # If you try to close the connection with
                          # "self._connection = None" in self.close(), after
                          # the connection was already closed, you have to make
                          # sure to not try to lookup '_connection' in
                          # self._status, which would fail, because the db is
                          # already closed, but self.is_open still returns a
                          # True.
                and name != '_connection'
                and name in self._status):
            self._status[name] = value
        else:
            super().__setattr__(name, value)

    def __str__(self):
        if self.filename is None:
            return self.__repr__()
        else:
            return ''.join([''.join(["'", self.filename, "':"]).ljust(20),
                            '<', self.__class__.__module__,
                            '.',
                            self.__class__.__name__, '>'])


def create_experiment(**kwargs):
    """
    Create a new Experiment or get the last created one.

    Parameters
    ----------
    **kwargs
        Parameters for the initialization of an instance of Experiment().

    Returns
    -------
    Experiment
        The new or the last created Experiment.
    """
    return Experiment.created_last or Experiment(**kwargs)


def open_experiment(filename=None, directory=None, **kwargs):
    """
    Create a new Experiment or get the last created one.
    Optionally, open an experiment file `filename'.

    Parameters
    ----------
    filename : str, optional
        The path to the file of the Experiment.
    directory : str, optional
        If `filename` alone is not sufficient to describe the location of the
        experiment file, you can provide a `directory`, where the path of
        `filename` is located in.
    **kwargs
        Parameters for the initialization of an instance of Experiment().

    Returns
    -------
    Experiment
        The Experiment with the opened file.
    """
    experiment = create_experiment(**kwargs)
    experiment.open(filename=filename, directory=directory)
    return experiment


def close_experiment():
    """
    Close the last created Experiment.
    """
    experiment = Experiment.created_last
    if experiment is not None:
        return experiment.close()
    return True


def save_experiment(**kwargs):
    """
    Save the last created Experiment.
    """
    experiment = Experiment.created_last
    if experiment is not None:
        return experiment.save(**kwargs)
    return True


def create_calibration(calibration_type=None, verbose=False, **kwargs):
    """
    Create a Calibration object.

    Parameters
    ----------
    calibration_type : str
        Select the type of the Calibration. Possible types can be found in the
        default `cfgfile` `pyoti/etc/calibration.cfg`
    verbose : bool
        Be verbose about the determination of the calibration_type upon reading
        the paramaters from the `cfgfile`.
    **kwargs
        Parameters will be passed to the function
        `calibration.create_calibration()`

    Returns
    -------
    Calibration
        A Calibration object.
    """
    calibration = cb.create_calibration(source_type=calibration_type,
                                        verbose=verbose, **kwargs)
    print('Created calibration with: ')
    print('  object: ' + str(calibration.calibsource))
    print('  name: ' + calibration.calibsource.name)
    print('  beta: ' + str(calibration.beta()))
    print('  kappa: ' + str(calibration.kappa()))
    print('If this seems suspicious to you, recheck the parameters you '
          'provided!')
    return calibration
