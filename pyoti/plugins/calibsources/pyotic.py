# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:32:58 2016

@author: Tobias Jachowski
"""
import numpy as np
import pint

from pyoti import config as cfg
from pyoti.calibration.calibsource import CalibrationSource
from pyoti import helpers as hp


class PyOTICfSource(CalibrationSource):
    def __init__(self, filename=None, focalshift=0.785, name=None, **kwargs):
        if filename is None:
            raise TypeError("PyOTICfSource missing the required positional "
                            "argument 'filename'.")
        self.filename = filename

        value, factor = read_pyotic_hc_result(filename)

        self.corrfactor = value['correction_factor']
        self.dsurf = value['apparent_surface_distance'] \
            * factor['apparent_distance']
        self.beta = np.array([value['beta_x'],
                              value['beta_y'],
                              value['beta_z']]) * factor['beta']
        self.kappa = np.array([value['kappa_x'],
                               value['kappa_y'],
                               value['kappa_z']]) * factor['kappa']
        self.mbeta = np.array([value['mbeta_x'],
                               value['mbeta_y'],
                               value['mbeta_z']]) * factor['mbeta']
        self.mkappa = np.array([value['mkappa_x'],
                                value['mkappa_y'],
                                value['mkappa_z']]) * factor['mkappa']
        self.radiusspec = value['radius'] * factor['radius']
        self.focalshift = focalshift
        self.name = name or 'OT Investigator and Calibrator Height ' \
                            'Calibration file originally loaded from \n    ' \
                            '%s' % (self.filename)


def read_pyotic_hc_result(filename):
    pars = cfg.read_cfg_file(filename)

    ureg = pint.UnitRegistry()

    results = pars['RESULTS']
    units = pars['UNITS']

    values = {}
    for key, value in results.items():
        try:
            values[key] = float(value)
        except:
            values[key] = value

    factors = {}
    factors['radius'] = ureg(units['radius']).to('um').magnitude
    factors['apparent_distance'] \
        = ureg(units['apparent_distance']).to('um').magnitude
    factors['beta'] = ureg(units['beta']).to('nm/mV').magnitude
    factors['mbeta'] = ureg(units['mbeta']).to('nm/(mV*um)').magnitude
    factors['kappa'] = ureg(units['kappa']).to('pN/nm').magnitude
    factors['mkappa'] = ureg(units['mkappa']).to('pN/(nm*um)').magnitude

    return values, factors


class PyOTICSource(CalibrationSource):
    def __init__(self, heightcalibration, focalshift=0.785, name=None,
                 **kwargs):
        self.hc = heightcalibration

        self.name = name or 'OT Investigator and Calibrator Height Calibration'

    @property
    def corrfactor(self):
        return self.hc.height_fit_results.corr
        # value['correction_factor']

    @property
    def dsurf(self):
        return self.hc.height_fit_results.h0
        # value['apparent_surface_distance']
        # * ureg(units['apparent_distance']).to('um').magnitude

    @property
    def beta(self):
        hfr = self.hc.height_fit_results
        return np.array([hfr.beta_x, hfr.beta_y, hfr.beta_z])
        # * ureg(units['beta']).to('nm/mV').magnitude

    @property
    def kappa(self):
        hfr = self.hc.height_fit_results
        return np.array([hfr.kappa_x, hfr.kappa_y, hfr.kappa_z])
        # * ureg(units['kappa']).to('pN/nm').magnitude

    @property
    def mbeta(self):
        hfr = self.hc.height_fit_results
        return np.array([hfr.mbeta_x, hfr.mbeta_y, hfr.mbeta_z])
        # * ureg(units['mbeta']).to('nm/(mV*um)').magnitude

    @property
    def mkappa(self):
        hfr = self.hc.height_fit_results
        return np.array([hfr.mkappa_x, hfr.mkappa_y, hfr.mkappa_z])
        # * ureg(units['mkappa']).to('pN/(nm*um)').magnitude

    @property
    def radiusspec(self):
        return self.hc.height_fit_results.radius
        # value['radius']
        # * ureg(units['radius']).to('um').magnitude

    @property
    def focalshift(self):
        return self.hc.focal_shift


class PyOTICSingleSource(CalibrationSource):
    def __init__(self, filename=None, number=1, active=True, radiusspec=None,
                 corrfactor=1.0, focalshift=None, mbeta=None, mkappa=None,
                 dsurf=0.0, name=None, **kwargs):
        if filename is None:
            raise TypeError("PyOTICSingle missing the required positional "
                            "argument 'filename'.")
        self.filename = filename
        self.number = 1
        self.active = active

        # read in the calibration data and parameters
        data, pars = read_pyotic_single_calibration(filename, number=number)

        # calculate conversion factors
        ureg = pint.UnitRegistry()
        factor = {}
        factor['radius'] = ureg(pars['radius_unit']).to('um').magnitude
        factor['beta'] \
            = ureg(pars['displacement_sensitivity_unit']).to('nm/mV').magnitude
        factor['kappa'] \
            = ureg(pars['trapstiffness_unit']).to('pN/nm').magnitude

        # extract the calibration data
        self.excited_axis = int(data['ex_axis'])
        beta_x_pc = data['beta_x_pc'] * factor['beta']
        beta_y_pc = data['beta_y_pc'] * factor['beta']
        beta_z_pc = data['beta_z_pc'] * factor['beta']
        self.beta_pc = np.array([beta_x_pc, beta_y_pc, beta_z_pc])
        beta_x_ac = data['beta_x_ac'] * factor['beta']
        beta_y_ac = data['beta_y_ac'] * factor['beta']
        beta_z_ac = data['beta_z_ac'] * factor['beta']
        self.beta_ac = np.array([beta_x_ac, beta_y_ac, beta_z_ac])
        kappa_x_pc = data['kappa_x_pc'] * factor['kappa']
        kappa_y_pc = data['kappa_y_pc'] * factor['kappa']
        kappa_z_pc = data['kappa_z_pc'] * factor['kappa']
        self.kappa_pc = np.array([kappa_x_pc, kappa_y_pc, kappa_z_pc])
        kappa_x_ac = data['kappa_x_ac'] * factor['kappa']
        kappa_y_ac = data['kappa_y_ac'] * factor['kappa']
        kappa_z_ac = data['kappa_z_ac'] * factor['kappa']
        self.kappa_ac = np.array([kappa_x_ac, kappa_y_ac, kappa_z_ac])
        if mbeta is None:
            mbeta = [0.0 for i in range(self.beta.size)]
        self.mbeta = np.array(hp.listify(mbeta))
        if mkappa is None:
            mkappa = [0.0 for i in range(self.kappa.size)]
        self.mkappa = np.array(hp.listify(mkappa))
        self.dsurf = dsurf
        self.radiusspec = radiusspec \
            or pars['specified_radius'] * factor['radius']
        self.corrfactor = corrfactor
        self.focalshift = focalshift or pars['focal_shift']
        self.name = name or 'OT Investigator and Calibrator single ' \
                            'calibration originally loaded from \n    %s' \
                            % (self.filename)

    @property
    def beta(self):
        if self.active:
            return self.beta_ac
        else:
            return self.beta_pc

    @beta.setter
    def beta(self, beta):
        print('Setting beta not supported. Instead, set the property '
              '`beta_pc` or `beta_ac`.')

    @property
    def kappa(self):
        if self.active:
            return self.kappa_ac
        else:
            return self.beta_pc

    @kappa.setter
    def kappa(self, kappa):
        print('Setting kappa not supported. Instead, set the property '
              '`kappa_pc` or `kappa_ac`.')


def read_pyotic_single_calibration(hc_datafile, hc_parafile=None, number=1):
    # Autodetermine the name of the parameterfile
    hc_parafile = hc_parafile or hc_datafile.replace('_data.dat',
                                                     '_parameters.txt')

    # read in the parameters of the calibration
    pars = cfg.read_cfg_file(hc_parafile)
    # convert strings into floats
    parameters = {}
    for key, value in pars['DEFAULT'].items():
        try:
            parameters[key] = float(value)
        except:
            parameters[key] = value

    # read in the data of the first calibration
    with open(hc_datafile) as f:
        lines = f.readlines()
        header = lines[0].split()
        calibration = lines[number].split()
    # create a dict with all the data of the calibration
    data = dict(zip(header, [float(i) for i in calibration]))

    return data, parameters
