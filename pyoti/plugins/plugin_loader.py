#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 18:13:42 2017

@author: tobiasj
"""

import importlib
import sys
import os
import appdirs

dirs = appdirs.AppDirs('pyoti')

__pyoti_plugins_dir__= os.path.dirname(globals()['__file__'])
__site_plugins_dir__ =  os.path.join(dirs.site_data_dir, 'plugins')
__user_plugins_dir__ = os.path.join(dirs.user_data_dir, 'plugins')
__plugins_dirs__ = [__pyoti_plugins_dir__, __site_plugins_dir__,
                    __user_plugins_dir__]
__plugin_kinds__ = ['datasources', 'modifications', 'calibsources']


def load_module(module_path, module_name):
    # Load module from source file [1], [2], and [3].
    # Set the __name__ attribute of the module and the __module__ attribute of
    # all classes to 'module_name'.
    # The names of the classes of objects, consisting of their class and module
    # names, are stored in the ZODB database, permanently [4]!
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # Map old modules to new module names that have already been loaded [5].
    # Objects are unpickled by looking up the classes of the corresponding
    # module.
    sys.modules[module_name] = module

    # [1] https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    # [2] https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path#67692
    # [3] https://stackoverflow.com/questions/19009932/import-arbitrary-python-source-file-python-3-3#19011259
    # [4] http://www.zodb.org/en/latest/tutorial.html#introduction
    # [5] https://docs.python.org/3/library/sys.html


def load_modules():
    for plugins_dir in __plugins_dirs__:
        for kind in __plugin_kinds__:
            path = os.path.join(plugins_dir, kind)
            if os.path.isdir(path):
                for name in os.listdir(path):
                    if name[-3:] != '.py' or name == '__init__.py':
                        continue
                    mod_path = os.path.join(plugins_dir, kind, name)
                    mod_name = ''.join(['pyoti.plugins.',
                                        kind, '.',
                                        name[:-3]])
                    load_module(mod_path, mod_name)