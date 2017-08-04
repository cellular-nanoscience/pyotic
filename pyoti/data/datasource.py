# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 19:13:51 2016

@author: Tobias Jachowski
"""
import os
import persistent
from functools import wraps

from .. import helpers as hp


def if_needs_file(func):
    @wraps(func)
    def func_wrapper(self, *args, **kwargs):
        if self._needs_file:
            return func(self, *args, **kwargs)
        print("DataSource does not need files to get data!")
    return func_wrapper


class DataSource(persistent.Persistent):
    name = 'Generic data source'
    samplingrate = 1000  # default 1000 Hz
    _needs_file = False

    def __init__(self, filename, directory=None, root=None, **kwargs):
        if filename is not None:
            self._needs_file = True
            self._root = root
            filename, absdir, fullfilename = hp.file_and_dir(filename,
                                                             directory)
            self._filename = filename
            self._directory = os.path.relpath(absdir, self.rootdir)
            self._filename_orig = filename
            self._directory_orig = self._directory
            self._absdir_orig = absdir

    def as_array():
        """
        This function must return a 2D numpy.ndarray, containing the data.
        """
        return None

    @property
    @if_needs_file
    def rootdir(self):
        """
        The directory of the experiment file.
        """
        if self._root is None:
            print("The datasource is not stored in an experiment file!\n"
                  "Therefore, no directory is known!\n"
                  "Instead, the current working directory is used!")
            return os.path.realpath('.')
        return self._root._v_absdir

    @property
    @if_needs_file
    def directory(self):
        """
        The directory of the data file relative to the experiment file.
        """
        return self._directory

    @property
    @if_needs_file
    def absdir(self):
        """
        The absolute directory of the data file.
        """
        return os.path.realpath(os.path.join(self.rootdir, self.directory))

    @property
    @if_needs_file
    def filename(self):
        """
        The filename of the data file.
        """
        return self._filename

    @property
    @if_needs_file
    def absfile(self):
        """
        The absolute path to the data file.
        """
        return os.path.join(self.absdir, self.filename)

    @property
    @if_needs_file
    def directory_orig(self):
        """
        The relative directory the data file was originally loaded from.
        """
        return self._directory_orig

    @property
    @if_needs_file
    def filename_orig(self):
        """
        The original name of the data file.
        """
        return self._filename_orig

    @property
    @if_needs_file
    def absdir_orig(self):
        """
        The absolute directory the data file was loaded from.
        """
        return self._absdir_orig

    @property
    @if_needs_file
    def absfile_orig(self):
        """
        The absolute path to the data file it was originally loaded from.
        """
        return os.path.join(self.absdir_orig, self.filename_orig)
