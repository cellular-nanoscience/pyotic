# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:41:17 2016

@author: Tobias Jachowski
"""
import inspect
import numbers

from pyoti.data.datasource import DataSource
from pyoti.picklable import unboundfunction


class GenericDataFile(DataSource):
    def __init__(self, load_data, filename, directory=None, samplingrate=1.0,
                 **kwargs):
        """
        load_data : function
        filename : str
        directory : str
        samplingrate : float
        **kwargs
        """
        super().__init__(filename=filename, directory=directory, **kwargs)

        self.load_data = unboundfunction(load_data)

        if isinstance(samplingrate, numbers.Number):
            self.samplingrate = samplingrate
        else:
            samplingrate_args = {}
            for par in inspect.getargspec(samplingrate)[0]:
                if par in kwargs:  # par can be anything, except load_data,
                                   # filename, directory or samplingrate
                    samplingrate_args[par] = kwargs.get(par)
                if par == 'filename':  # automatically use filename
                    samplingrate_args['filename'] = self.absfile
            self.samplingrate = samplingrate(**samplingrate_args)

        self.load_data_args = {}
        for par in inspect.getargspec(load_data)[0]:
            if par in kwargs:
                self.load_data_args[par] = kwargs.get(par)

        self.name = ("Generic data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (self.absfile_orig,
                                                  self.samplingrate)

    def as_array(self):
        filename = self.absfile
        data = self.load_data(filename, **self.load_data_args)
        return data


class GenericData(DataSource):
    def __init__(self, load_data, samplingrate=1.0, **kwargs):
        """
        load_data : function
        samplingrate : float
        """
        self.load_data = unboundfunction(load_data)

        if isinstance(samplingrate, numbers.Number):
            self.samplingrate = samplingrate
        else:
            samplingrate_args = {}
            for par in inspect.getargspec(samplingrate)[0]:
                if par in kwargs:
                    samplingrate_args[par] = kwargs.pop(par)
            self.samplingrate = samplingrate(**samplingrate_args)

        self.fun_args = {}
        for par in inspect.getargspec(load_data)[0]:
            if par in kwargs:
                self.fun_args[par] = kwargs.pop(par)

        self.name = ("Generic data with \n"
                     "    samplingrate %s Hz") % (self.samplingrate)

    def as_array(self):
        data = self.load_data(**self.fun_args)
        return data
