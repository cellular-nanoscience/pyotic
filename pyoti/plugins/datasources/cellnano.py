# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:37:58 2016

@author: Tobias Jachowski
"""
import numpy as np
import os

from pyoti.data import labview as lv
from pyoti.data.datasource import DataSource


class CNLabViewBinData(DataSource):
    def __init__(self, filename, directory=None, parafile=None, datext='.bin',
                 parext='_para.dat', **kwargs):
        """
        parext : str, optional
            The extension of the parameter file ('_para.dat', default)
        """
        super().__init__(filename=filename, directory=directory, **kwargs)

        if parafile is None:
            parafile = self.absfile.replace(datext, parext)
        self._parafile_orig = parafile

        self.samplingrate = get_samplingrate_from_parameter_file(parafile)

        self.name = ("LabVIEW bin data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (self.absfile_orig,
                                                  self.samplingrate)

    def as_array(self):
        filename = self.absfile
        data = lv.read_labview_bin_data(filename)
        return data

    @property
    def parafile_orig(self):
        return self._parafile_orig


class CNLabViewTxtData(DataSource):
    def __init__(self, filename, directory=None, samplingrate=1000.0,
                 **kwargs):
        """
        parext : str, optional
            The extension of the parameter file ('_para.dat', default)
        """
        super().__init__(filename=filename, directory=directory, **kwargs)
        self.samplingrate = samplingrate

        self.name = ("LabVIEW txt data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (self.absfile_orig,
                                                  self.samplingrate)

    def as_array(self):
        filename = self.absfile
        data = np.loadtxt(filename, skiprows=5)
        return data


class CNLabViewQPDviData(CNLabViewTxtData):
    def __init__(self, filename, directory=None, **kwargs):
        super().__init__(filename, directory=directory, **kwargs)
        self.samplingrate = np.genfromtxt(self.absfile, skip_header=2,
                                          max_rows=1)[0]

        self.name = ("LabVIEW QPD.vi data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (os.path.join(self.directory,
                                                               self.filename),
                                                  self.samplingrate)


def get_samplingrate_from_parameter_file(parafile):
    para = np.loadtxt(parafile, comments='%', delimiter='\t')
    scanrate = para[0]  # float
    decimating = para[2]  # float
    samplingrate = scanrate / decimating  # float
    return samplingrate
