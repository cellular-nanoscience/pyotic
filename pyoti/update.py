# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 18:37:33 2016

@author: Tobias Jachowski
"""
import collections
import os

from . import helpers as hp
from . import version
from .picklable import Attributes


def comp(version):
    return tuple(map(int, (version.split('.'))))


def update_db(experiment):
    """
    Update the format of the file the Experiment is saved in, if the structure
    of the database has changed.
    """
    if not experiment.is_open:
        print("Experiment is closed!")
        return

    old_version = experiment.version
    new_version = version()

    _update_db(experiment, old_version, new_version)


def _update_db(experiment, old_version, new_version):
    if comp(old_version) < comp(new_version):
        print("Updating database from %s to %s ..." % (old_version,
                                                       new_version))

        old_version = comp(old_version)

        # Load analyzer module from tethered_bead package (<= 0.2.3)
        if old_version <= comp('0.2.3'):
            from tethered_bead import analyzer
            # prevent modifications from trying to load data upon loading
            analyzer.Modification.__setstate__ \
                = analyzer.GraphMember.__setstate__

        # Load investigator package (<= 0.5.2)
        if old_version <= comp('0.5.2'):
            import investigator
            # prevent modifications from trying to load data upon loading
            investigator.modification.Modification.__setstate__ \
                = investigator.graph.GraphMember.__setstate__

            # Update all objects to 0.5.1+ (0.5.2), where the module name
            # was changed from investigator to pyoti. If you import the
            # package from the new path only and try to open a database
            # stored with the old package name, you will get broken
            # objects.
            # This means, you have to import the package from the old AND
            # the new path at the same time and then convert all instances
            # that have their classes defined in the investigator package:
            import importlib

            # investigator and pyoti need to be imported ahead
            import pyoti

            def change_object_module(obj, old='investigator', new='pyoti',
                                     new_class_name=None):
                class_name = obj.__class__.__name__
                old_module_name = obj.__module__
                new_module_name = old_module_name.replace(old, new, 1)
                new_module = importlib.import_module(new_module_name)
                if new_class_name is not None:
                    class_name = new_class_name
                obj.__class__ = getattr(new_module, class_name)
                if hasattr(obj, '_p_changed'):
                    obj._p_changed = True

        # Update to 0.2.2
        if old_version < comp('0.2.2'):
            # Update to 0.2.3
            for record in experiment._graphroot.children:
                calibsource = record.calibration.calibsource
                record._traces = ['psdX', 'psdY', 'psdZ',
                                  'positionX', 'positionY', 'positionZ',
                                  'mirrorX', 'laser']
                if hasattr(calibsource, 'C'):
                    calibsource.beta = calibsource.C
                    del calibsource.C
                if hasattr(calibsource, 'mC'):
                    calibsource.mbeta = calibsource.mC
                    del calibsource.mC
                if hasattr(calibsource, 'k'):
                    calibsource.kappa = calibsource.k
                    del calibsource.k
                if hasattr(calibsource, 'mk'):
                    calibsource.mkappa = calibsource.mk
                    del calibsource.mk
                if hasattr(calibsource, 'hsurf'):
                    calibsource.dsurf = calibsource.hsurf
                    del calibsource.hsurf
                # if hasattr(record, '_inversion'):
                    # positionZ was inverted, it now is proportional to the
                    # distance of the bead to the surface and not the
                    # position of the stage anymore.
                    # record._inversion[5] = -1
            for modification in experiment._graphroot.members(
                    instance_class=analyzer.Impact):
                if 'stageZ' in modification._traces_apply:
                    traces = modification._traces_apply
                    i = traces.index('stageZ')
                    traces.remove('stageZ')
                    traces.insert(i, 'positionZ')
                if hasattr(modification, 'wa'):
                    del modification.wa
                if hasattr(modification, 'leftpull'):
                    del modification.leftpull
                if hasattr(modification, 'rightpull'):
                    del modification.rightpull
            for modification in experiment._graphroot.members(
                    instance_class=analyzer.Rotation):
                if hasattr(modification, 'wa'):
                    del modification.wa
                if hasattr(modification, 'stageZ'):
                    del modification.stageZ
            for modification in experiment._graphroot.members(
                    instance_class=analyzer.Attachment):
                if hasattr(modification, 'wa'):
                    del modification.wa
                if hasattr(modification, 'rightstage'):
                    del modification.rightstage
                if hasattr(modification, 'leftstage'):
                    del modification.leftstage
                if hasattr(modification, 'pullleft'):
                    del modification.pullleft
                if hasattr(modification, 'pullright'):
                    del modification.pullright

        # Update to 0.1.3
        if old_version < comp('0.1.3'):
            # Update to version 0.1.3
            for record in experiment._graphroot.children:
                if hasattr(record, 'datasource') \
                        and not hasattr(record.datasource, '_datadir_orig'):
                    try:
                        if hasattr(record.datasource, 'datadir'):
                            datasource = record.datasource
                            datasource._datadir_orig \
                                = record.datasource.datadir
                        if isinstance(record.datasource,
                                      analyzer.ASWADBinData):
                            datasource.name = (
                                "ASWAD LabVIEW bin data originally from \n"
                                "    %s%s%s with \n"
                                "    samplingrate %s Hz") \
                                % (datasource.datadir_orig,
                                   os.sep,
                                   datasource.datafile,
                                   datasource.samplingrate)
                    except:
                        pass

        # Update to 0.1.4
        if old_version < comp('0.1.4'):
            # Update to version 0.2.0
            for modification in experiment.modifications():
                if isinstance(modification, analyzer.Rotation):
                    del modification.stageZ
                    del modification.calibration

        # Update to 0.2.3 database format
        if old_version <= comp('0.2.3'):
            # Retrieve all objects in database ...
            change_object_module(experiment._graphroot,
                                 old='tethered_bead.analyzer',
                                 new='pyoti.graph')
            for obj in experiment._graphroot.members():
                if isinstance(obj, analyzer.Record):
                    # Record
                    change_object_module(obj,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.region.record')
                    # Node of calibration
                    change_object_module(obj.calibration._node,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.graph')
                    # Calibration
                    change_object_module(obj.calibration,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.calibration.calibration')
                    # CalibSource
                    calibsource = obj.calibration.calibsource
                    if isinstance(calibsource, analyzer.CalibrationMATLAB):
                        calibsource.filename = calibsource.calibfile
                        del calibsource.calibfile
                        change_object_module(
                            calibsource,
                            old='tethered_bead.analyzer',
                            new='pyoti.calibration.calibsources.cnmatlab',
                            new_class_name='CNMSource')
                    elif isinstance(calibsource, analyzer.CalibrationTA):
                        calibsource.filename = calibsource.tacfile
                        del calibsource.tacfile
                        change_object_module(
                            calibsource,
                            old='tethered_bead.analyzer',
                            new='pyoti.calibration.calibsources.cnmatlab',
                            new_class_name='CNMSource')
                    else:
                        change_object_module(
                            calibsource,
                            old='tethered_bead.analyzer',
                            new='pyoti.calibration.calibsource')
                    # DataSource
                    datasource = obj.datasource
                    if isinstance(datasource, analyzer.ASWADBinData):
                        change_object_module(
                            datasource,
                            old='tethered_bead.analyzer',
                            new='pyoti.data.datasources.cnlabview',
                            new_class_name='CNLabViewBinData')
                        # Convert ASWADBinData to CNLabViewBinData
                        filename, absdir, fullfilename \
                            = hp.file_and_dir(datasource.datafile,
                                              experiment._graphroot._v_absdir)
                        datasource._filename = filename
                        datasource._directory \
                            = os.path.relpath(absdir,
                                              experiment._graphroot._v_absdir)
                        datasource._filename_orig = filename
                        datasource._directory_orig \
                            = datasource._datadir_orig
                        _fn, datasource._absdir_orig, _ffn \
                            = hp.file_and_dir(datasource.datafile,
                                              datasource._datadir_orig)
                        del datasource.datadir
                        del datasource.datafile
                        del datasource._datadir_orig
                # Node
                change_object_module(obj._node,
                                     old='tethered_bead.analyzer',
                                     new='pyoti.graph')
                # MultiRegion
                if isinstance(obj, analyzer.MultiRegion):
                    change_object_module(obj,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.region.view')
                # View
                if isinstance(obj, analyzer.View):
                    change_object_module(obj,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.region.view')
                # InteractiveAttributes of Modification
                if isinstance(obj, analyzer.Modification):
                    change_object_module(obj.iattributes,
                                         old='tethered_bead.analyzer',
                                         new='pyoti.picklable')
                # Offset modifiation
                if isinstance(obj, analyzer.Offset):
                    change_object_module(
                        obj,
                        old='tethered_bead.analyzer',
                        new='pyoti.modification.modifications.offset')
                # Impact modification
                if isinstance(obj, analyzer.Impact):
                    change_object_module(
                        obj,
                        old='tethered_bead.analyzer',
                        new='pyoti.modification.modifications.impact')
                # Beadscan modification
                if isinstance(obj, analyzer.Beadscan):
                    change_object_module(
                        obj,
                        old='tethered_bead.analyzer',
                        new='pyoti.modification.modifications.beadscan')
                # Attachment modification
                if isinstance(obj, analyzer.Attachment):
                    change_object_module(
                        obj,
                        old='tethered_bead.analyzer',
                        new='pyoti.modification.modifications.attachment')
                # Rotation modification
                if isinstance(obj, analyzer.Rotation):
                    change_object_module(
                        obj,
                        old='tethered_bead.analyzer',
                        new='pyoti.modification.modifications.rotation')
                # Modification boundmethod and evaluate button
                if isinstance(obj,
                              pyoti.modification.modification.Modification):
                    # picklable_boundmethod -> pyoti.picklable.boundmethod
                    for functions in \
                            obj.iattributes._callback_functions.values():
                        for fun in functions:
                            change_object_module(fun,
                                                 old='tethered_bead.analyzer',
                                                 new='pyoti.picklable',
                                                 new_class_name='boundmethod')
                    # Add evaluate button
                    obj.add_iattribute('evaluate', description='Evaluate',
                                       unset_automatic=False,
                                       set_changed=False,
                                       callback_functions=[obj.evaluate])

        # Switch off inversion of 'positionZ' before version 0.2.3
        if old_version <= comp('0.2.3'):
            for record in experiment._graphroot.children:
                record.set_inversion('positionZ', False)

        # Update to 0.5.0
        if old_version < comp('0.5.0'):
            # The database structure for iattributes has changed.
            # Modifications are the only objects using iattributes.
            # Update all iattributes by triggering their unpickling and
            # informing ZODB of change
            # The update itself is done in the method:
            # picklable.InteractiveAttributes.__setstate__()
            for modification in experiment.modifications():
                modification.iattributes._p_changed = True

        # Update to 0.5.2 database format
        if old_version <= comp('0.5.2'):
            if experiment.cached_region is None:
                experiment.cached_region = []
            if not isinstance(experiment.cached_region, collections.Iterable):
                experiment.cached_region = [experiment.cached_region]

            # Retrieve all objects in database ...
            for obj in experiment._graphroot.members():
                if isinstance(obj, pyoti.region.Record) \
                        or isinstance(obj, investigator.region.Record):
                    # Node of calibration
                    change_object_module(obj.calibration._node)
                    # Calibration
                    change_object_module(obj.calibration)
                    # CalibSource
                    change_object_module(obj.calibration.calibsource)
                    # DataSource
                    change_object_module(obj.datasource)
                if isinstance(obj, pyoti.modification.Modification) \
                    or isinstance(obj,
                                  investigator.modification.Modification):
                    change_object_module(obj.iattributes)
                # Node
                change_object_module(obj._node)
                # Cargo (Record, MultiRegion, View, Modification)
                change_object_module(obj)

            # change_object_module(experiment._dbroot['status'])
            # experiment._dbroot._p_changed = True
            new_status = Attributes()
            old_status = experiment._dbroot['status']
            new_status.cached_region = old_status.cached_region
            if hasattr(old_status, 'last_added_region'):
                new_status.last_added_region = old_status.last_added_region
            else:
                new_status.last_added_region = None
            if hasattr(old_status, 'last_adjusted_view'):
                new_status.last_adjusted_view = old_status.last_adjusted_view
            else:
                new_status.last_adjusted_view = None
            experiment._dbroot['status'] = new_status

        # Update to 0.5.3
        if old_version < comp('0.5.3'):
            for child in experiment._graphroot.children:
                if (hasattr(child, 'datasource')
                        and hasattr(child.datasource, 'directory')):
                    try:
                        child.datasource._load_file = True
                        child.datasource._root = experiment._graphroot
                        child.datasource._filename \
                            = child.datasource._filename_orig
                        child.datasource._directory = '.'
                        # Unfortunately, we can't convert the relative
                        # directory_orig to an absolute path, because we don't
                        # know, whether the file is opened from the same
                        # location, during the update.
                        child.datasource._absdir_orig \
                            = child.datasource._directory_orig
                        del child.datasource.__dict__['directory']
                        del child.datasource.__dict__['filename']
                        # The last one could fail, but this doesn't matter.
                        del child.datasource.__dict__['parafile']
                    except:
                        pass

        # Update to 0.5.4
        if old_version < comp('0.5.4'):
            for region in experiment.regions():
                if hasattr(region, '_data_cached'):
                    del region.__dict__['_data_cached']
                    region._p_changed = True

        # Update to 0.6.0
        if old_version < comp('0.6.0'):
            for view in experiment.views():
                sticky = True
                if hasattr(view, 'autonomous'):
                    sticky = not view.autonomous
                    del view.__dict__['autonomous']
                view.sticky_start = sticky
                view.sticky_stop = sticky

        # Update to 0.6.1
        if old_version < comp('0.6.1'):
            for mod in experiment.modifications():
                try:
                    if hasattr(mod, 'ifigure'):
                        del mod.__dict__['ifigure']
                        mod._p_changed = True
                except:
                    pass
                try:
                    if hasattr(mod, 'rfigure'):
                        del mod.__dict__['rfigure']
                        mod._p_changed = True
                except:
                    pass

        # Update to 0.6.2
        if old_version < comp('0.6.2'):
            for child in experiment._graphroot.children:
                if hasattr(child, 'datasource'):
                    ds = child.datasource
                    # force datasource to be unpickled
                    if hasattr(ds, '_load_file') \
                            and '_load_file' in ds.__dict__:
                        ds._needs_file = ds.__dict__['_load_file']
                        del ds.__dict__['_load_file']

        # Update to 0.8.0
        if old_version < comp('0.8.0'):
            # Convert all GModifications to corresponding new Modifications
            # GAttachment -> Attachment
            # GBaseline-> Baseline
            # GBeadscan -> Beadscan
            #  Beadscan.spline -> Beadscan.model
            #  delete beadscan.bin_means
            # GImpact -> Touchdown
            #  pyoti.modification.modifications.impact
            #      -> pyoti.modification.modifications.touchdown
            #  Impact -> Touchdown
            #  iattributes.impact -> iattributes.touchdown
            #  iattributes.set_impact_calib -> iattributes.set_touchdown_calib
            #   delete impact.bin_means
            # GRotation -> Rotation
            # for all iattributes datapoints_used -> datapoints
            from .modification import \
                Attachment, Baseline, Beadscan, Touchdown, Rotation

            MODIFICATION_CLASS = {
                'GAttachment': Attachment,
                'GBaseline': Baseline,
                'GBeadscan': Beadscan,
                'GImpact': Touchdown,
                'GRotation': Rotation
            }

            def create_new_mod(old_mod, new_modclass):
                view_based = old_mod.view_based
                view_apply = old_mod.view_apply
                name = old_mod.name
                new_mod = new_modclass(view_based=view_based,
                                       view_apply=view_apply,
                                       name=name, group=name,
                                       db_update=True)
                new_mod._traces_apply = old_mod.traces_apply
                if hasattr(old_mod.iattributes, 'automatic'):
                    new_mod.iattributes.automatic \
                        = old_mod.iattributes.automatic
                return new_mod

            def set_iattribute_value(old_mod, new_mod, key, new_key=None):
                new_key = new_key or key
                value = old_mod.iattributes[key]
                new_mod.iattributes.set_value(new_key, value, callback=False)

            # create new mods
            for old_mod in experiment.modifications():
                old_class_name = old_mod.__class__.__name__
                try:
                    new_class = MODIFICATION_CLASS[old_class_name]
                except:
                    continue
                new_mod = create_new_mod(old_mod, new_class)

                # transfer individual attributes
                if old_class_name == 'GBaseline':
                    new_mod._baseline_idx = old_mod._baseline_idx.copy()
                    new_mod._model = old_mod.model.copy()
                if old_class_name == 'GBeadscan':
                    new_mod._model = old_mod.spline.copy()
                    new_mod.minpsdZ = old_mod.minpsdZ
                    new_mod.maxpsdZ = old_mod.maxpsdZ
                if old_class_name == 'GImpact':
                    new_mod.left_upper = old_mod.left_upper
                    new_mod.right_upper = old_mod.right_upper
                    new_mod.left = old_mod.left
                    new_mod.right = old_mod.right
                if old_class_name == 'GRotation':
                    new_mod.rotation_method = old_mod.rotation_method

                # transfer iattributes values
                for key in old_mod.iattributes._widgets:
                    new_key = key
                    if key == 'datapoints_used':
                        new_key = 'datapoints'
                    if key == 'impact':
                        new_key = 'touchdown'
                    if key == 'set_impact_calib':
                        new_key = 'set_touchdown_calib'
                    set_iattribute_value(old_mod, new_mod, key, new_key)

                # remove old mod
                old_mod.view_based.remove_child(old_mod)
                # Second, delete relation to child
                old_mod.view_apply.remove_parent(old_mod)

        # Update to 0.9.0
        if old_version < comp('0.9.0'):
            # remove '_mod' from the group attribute. Since 0.8.1 name and
            # group can be independently set.
            for member in experiment.members():
                member.group = member.group.replace('_mod', '')
                if hasattr(member, 'iattributes'):
                    widget = member.iattributes._widgets['active']
                    widget.description = widget.description.replace(
                        "modification ", "")
                    if "datapoints" in member.iattributes._widgets:
                        widget = member.iattributes._widgets['datapoints']
                        widget.description = widget.description.replace(
                            "used ", "")

        # Update to 0.11.0
        if old_version < comp('0.11.0'):
            # Map old modules to new module names (new modules have already
            # been loaded by plugin.plugin_loader.load_modules())
            import sys
            import pyoti
            old_mods = ['pyoti.modification.modifications.attachment',
                        'pyoti.modification.modifications.baseline',
                        'pyoti.modification.modifications.beadscan',
                        'pyoti.modification.modifications.offset',
                        'pyoti.modification.modifications.rotation',
                        'pyoti.modification.modifications.touchdown',
                        'pyoti.modification.modifications.generic',
                        'pyoti.data.datasources.cnlabview',
                        'pyoti.data.datasources.generic',
                        'pyoti.calibration.calibsources.cellnano',
                        'pyoti.calibration.calibsources.pyotic']
            new_mods = ['pyoti.plugins.modifications.attachment',
                        'pyoti.plugins.modifications.baseline',
                        'pyoti.plugins.modifications.beadscan',
                        'pyoti.plugins.modifications.offset',
                        'pyoti.plugins.modifications.rotation',
                        'pyoti.plugins.modifications.touchdown',
                        'pyoti.plugins.modifications.generic',
                        'pyoti.plugins.datasources.cellnano',
                        'pyoti.plugins.datasources.generic',
                        'pyoti.plugins.calibsources.cellnano',
                        'pyoti.plugins.calibsources.pyotic']
            for old, new in zip(old_mods, new_mods):
                sys.modules[old] = sys.modules[new]

        # Update to 0.10.1
        if old_version < comp('0.10.1'):
            # Rename CNMSource to CNMatlabSource
            # pyoti.calibration.calibsources.cnmatlab.CNMSource ->
            # pyoti.calibration.calibsources.cellnano.CNMatlabSource
            for record in experiment.records():
                old_calibsource = record.calibration.calibsource
                old_class_name = old_calibsource.__class__.__name__
                if old_class_name == 'CNMSource':
                    change_object_module(
                        old_calibsource,
                        old='pyoti.calibration.calibsources.cnmatlab',
                        new='pyoti.calibration.calibsources.cellnano',
                        new_class_name='CNMatlabSource')

        # Update to 0.11.0
        if old_version < comp('0.11.0'):
            # To update the __module__ attribute of objects, one must
            # trigger a change of objects directly referencing the objects,
            # whose __module__ has changed, not the objects themselves. This is
            # weird!
            for obj in experiment._graphroot.members():
                if isinstance(obj, pyoti.modification.Modification):
                    # Update modification (directly referenced by a node)
                    obj._node._p_changed = True
                if isinstance(obj, pyoti.region.Record):
                    # Update datasource (directly referenced by a record)
                    obj._p_changed = True
                    # Update calibsource (directly referency by a calibration)
                    obj.calibration._p_changed = True

        # Update to X.Y.Z
        # if old_version < comp('X.Y.Z'):
            # for record in experiment.records():
            # for modification in experiment._graphroot.members(
            #                                           instance_class=Impact):
            # for record in experiment._graphroot.children:
            # for modification in experiment.modifications():
            # experiment._dbroot['version'] = 'X.Y.Z'

        # If there were no database layout changes, automatically set db
        # version to the newest version
        if comp(new_version) > comp(experiment._dbroot['version']):
            experiment._dbroot['version'] = new_version
