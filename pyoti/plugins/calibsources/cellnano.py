# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:31:12 2016

@author: Tobias Jachowski
"""
import numpy as np

from pyoti.calibration.calibsource import CalibrationSource


class CNMatlabSource(CalibrationSource):
    def __init__(self, filename=None, **kwargs):
        if filename is None:
            raise TypeError("CNMatlabSource missing the required positional "
                            "argument 'filename'.")
        self.filename = filename

        calib = np.genfromtxt(self.filename, skip_header=6,
                              invalid_raise=False)
        self.corrfactor = calib[0]
        self.dsurf = calib[1]
        self.beta = calib[2:5]
        self.kappa = calib[5:8]
        self.mbeta = calib[8:11]
        self.mkappa = calib[11:14]
        self.radiusspec = calib[14]
        self.focalshift = calib[15]
        self.name = "MATLAB calibration file originally loaded from \n    %s" \
                    % (self.filename)


class CNParaSource(CalibrationSource):
    def __init__(self, filename=None, beta=None, kappa=None, name=None,
                 **kwargs):
        if filename is None:
            raise TypeError("CNParaSource missing the required positional "
                            "argument 'filename'.")
        self.filename = filename
        para = np.loadtxt(self.filename, comments='%', delimiter='\t')
        beta = beta or para[3:6]
        kappa = kappa or para[6:9]
        name = name or 'Cellular Nanoscience parameter file originally ' \
                       'loaded from \n    %s' % (self.filename)
        super().__init__(beta=beta, kappa=kappa, name=name, **kwargs)
