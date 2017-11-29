# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
"""
Manage height-dependent power spectral density measurements and their fits.
"""
import pdb

from . import ureg

from . import name_constants as co

from .data_parser import read_std_data_file

from .focal_shift import get_focal_shift

from .physics import faxen_factor

from .plotting import add_plot_to_figure
from .plotting import col_dict
from .plotting import get_gridded_axis
from .plotting import get_residual_plot_axes
from .plotting import set_sym_y_labels

from .psd import PSDMeasurement
from .psd import gen_PSD_from_time_series
from .psd import ExpSetting
from .psd_fitting import is_outlier
from .psd_fitting import PSDFit

from .utilities import flatten_list
from .utilities import gen_fit_pars
from .utilities import u2str

from collections import OrderedDict

import configparser

from lmfit import minimize
from lmfit import report_fit

import matplotlib.pyplot as plt

from os import listdir

from os.path import join

from scipy import array
from scipy import exp
from scipy import inf
from scipy import logical_or
from scipy import pi
from scipy import shape
from scipy import sin
from scipy import zeros
from scipy import sqrt

import time

import warnings


def _lin_fun(x, x0=0, y0=0, slope=1):
    """
    Return the function value of a linear function y = (x - x0) * slope + y0.

    Note
    ----
    y0 is value at x0
    """
    return (x - x0) * slope + y0


def gen_height_fit_pars(names=None,
                        h0=0.0,  # in um
                        h0_min=-inf,  # in um
                        h0_max=inf,  # in um
                        corr=1.0,
                        corr_min=0.0,
                        corr_max=10.0,
                        focal_shift=0.8,
                        focal_shift_min=0.0,
                        focal_shift_max=2.0,
                        pars_named={},
                        add_pars=True,
                        **kwargs):
    """
    Create the lmfit.Parameters object for the height fit.

    The parameters beta, mbeta, kappa and mkappa will get a suffix '_name'
    according to one of the strings in the names list.

    Arguments
    ---------
    names : list(str)
        Names of the axes used, usually ['x', 'y', 'z']
    h0 : float
        Offset to the real surface height, in micrometers.
    hmin, hmax : float
        Lower and upper bound of the height offset.
    corr : float
        Correction factor depending on the fit method, either correcting for
        radius or viscosity deviations.
    corrmin, corrmax : float
        Lower and upper bound of the correction factor.
    beta_, mbeta_ : float
        Intecept and slope of the displacement senitivity depending on height
        in units nm/mV.
    kappa_, mkappa_ : float
        Intercept and slope of the trap stiffness in units pN/nm.

    Returns
    -------
    lmfit.Parameters

    Notes
    -----
    pars_named is a dictionary whose keys privide the names for parameters, the
    values are tuples of (initialization value, lower_bound, upper_bound)
    """
    init_kw = {'h0': h0,
               'h0_min': h0_min,
               'h0_max': h0_max,
               'corr': corr,
               'corr_min': corr_min,
               'corr_max': corr_max,
               'focal_shift': focal_shift,
               'focal_shift_min': focal_shift_min,
               'focal_shift_max': focal_shift_max
               }

    pn = {'beta': (0.5, 0, 40),
          'mbeta': (0.001, -1, 1),
          'kappa': (0.3, 0, 40),
          'mkappa': (0.001, -1, 1)}
    if add_pars:
        pars_named.update(pn)

    suffixes = ['', 'min', 'max']

    for name in names if names else []:
        for par_named, value in pars_named.items():
            for idx, suffix in enumerate(suffixes):
                if idx < 1:
                    key = '_'.join([par_named, name])
                else:
                    key = '_'.join([par_named, name, suffix])
                if key in kwargs:
                    init_kw[key] = kwargs[key]
                else:
                    init_kw[key] = value[idx]
    pars = gen_fit_pars(**init_kw)

    return pars


def _rel_drag_fit_fun(parameters,
                      heights,
                      radius,
                      fit_data=None,
                      errors=None,
                      lateral=True,
                      method='viscosity'
                      ):
    """
    Residual function for the relative drag.

    If fit_data and errors are None, the values for the parameters are
    calculated. If they are provided the resiuduals are returned.

    Arguments
    ---------
    parameters : lmfit.Parameters
        If None, default paremeters are generated.
        h0 - offset of the height data
        corr - correction factor to account for deviations in viscosity or
               correction factor to account for deviations from the specified
               bead radius, depending on if method set to viscosity or radius,
               respectively.
    heights : array(float)
        Height vector for the corresponding drag values in um.
    radius : float
        Radius of the trapped particle in um.
    focal_shift : float
        Focal shift of the trap.
    fit_data : array
        Relative drag values to fit to.
    errors : array
        The errors associated to the relative drag values.
    lateral : bool
        Direction of the movement. Lateral is parallel to the surface,
    method : str
        Determines if the correction factor 'corr' is associated with the
        'radius' or the 'viscosity'.

    Returns
    -------
    array

    See Also
    --------
    focal_shift
    """
    h = heights
    r = radius

    p = parameters.valuesdict()
    try:
        h0 = p['h0']
        corr = p['corr']
        fs = p['focal_shift']
    except:
        raise Exception('Fitting parameters are '
                        'missing in Parameters() object.')

    l = fs * (h - h0)

    if method is 'viscosity':
        rad = r
    else:
        rad = corr * r

    result = array(corr * faxen_factor(l, rad, lateral=lateral))

    if fit_data is None:
        return result

    if errors is None:
        return (fit_data - result)
    else:
        return (fit_data - result) / errors


def fit_rel_drag(heights,
                 rel_drag,
                 errors,
                 radius,
                 parameters=None,
                 lateral=True,
                 method='viscosity',
                 fit_focal_shift=False,
                 verbose=False,
                 **init_kwargs):
    """
    Calculates a least squares fit to the relative drag data.
    The fitting function used is 'rel_drag_fun()':
        rel_drag = drag / drag_stokes
                 = c_visc * c_rad * FL(c_rad*radius, height)

    Arguments
    ---------
    heights : array(float)
        Height vector for the corresponding drag values in micrometers.
    rel_drag : array(float)
        Relative drag values.
    errors : array(float)
        Error of the ralative drag values.
    focal_shift : float
        Factor that describes the physical shift of the trap above the surface.
    parameters : lmfit.Parameters
        If None, default paremeters are generated.
        h0 - offset of the height data
        corr - correction factor to account for deviations in viscosity or
               correction factor to account for deviations from the specified
               bead radius, depending on if method set to viscosity or radius,
               respectively.
    lateral : bool
        Specifies the Faxen's law correction of the drag. True referes to
        the axes parallel and False orthogonal to the surface.
    method : {'viscosity', 'radius'}
        whether correcting for viscosity or radius.
        N.B. One can use the method 'viscosity' still to correct
        for deviations of the radius, if the deviations are small and one can
        be sure that the viscosity did not deviate, e.g. due to local heating
        or other ingredients in the solution.
    verbose : bool
        Print some info to stdio.
    **init_kwargs : keyword arguments
        Keyword arguments defining the initial values of the fit, These are
        passed to gen_height_fit_pars()
        e.g. h0 = 0.0, corr = 1.0, h0min = -1.0, etc.

    Note
    ----
    The idea behind the drag_fit_method is to account for possible deviations
    from the 'known' viscosity or the radius, which would result in a false
    normalization of the drag to the stokes drag.
    I.e.:
        drag = drag_stokes' * faxen's law
             = 6 * pi * viscosity' * radius' * FL(radius', height)
             = 6 * pi * c_visc * viscosity * c_rad * radius
                 * FL(c_rad * radius, height)
             = drag_stokes * c_visc * c_rad * FL(c_rad * radius, height)
    where viscosity (of water) can be calculated if the temperature is known
    and radius is known from the bead manufacturer. Thus, one can write:
        rel_drag = drag / drag_stokes
                 = c_visc * c_rad * FL(c_rad*radius, height)

    The idea is now to set one of the two correction factors fixed to 1.0 and
    let the other factor being varied within the fitting algorithm. The method
    specifies which correction factor is allowed to vary.
    """
    r = radius

    if parameters is None:
        if 'corr' not in init_kwargs:
            init_kwargs['corr'] = min(rel_drag)
        pars = gen_height_fit_pars(**init_kwargs)
    else:
        pars = parameters

    if fit_focal_shift:
        pars['focal_shift'].vary = True
    else:
        pars['focal_shift'].vary = False

    if verbose:
        print('Fit to relative drag')

    minimizer = minimize(_rel_drag_fit_fun, pars,
                         args=(heights, r, rel_drag, errors, lateral, method),
                         nan_policy='propagate')
    if minimizer.success:
        if verbose:
            print('Least squares fit to relative drag succeeded.')
            report_fit(minimizer)
    else:
        if verbose:
            print('Least squares fit to relative drag did not succeed.')
            print('Minimizer message: "{0:s}"'.format(minimizer.message))
            print('lmdif message: "{0:s}"'.format(minimizer.lmdif_message))
    return minimizer


def _dissens_fit_fun(parameters,
                     height,  # in um
                     dissens_surf,  # in nm/mV
                     slope,  # in nm/(mV*um)
                     radius,  # in um
                     lateral=True,
                     method='viscosity',
                     fit_osci=False,
                     fun_kw={}  # ref_ind and wavelength (um)
                     ):
    """
    Return the value(s) of the theoretically acquired experimental
    displacment sensitivity, at a given height above the surface and
    neglecting the presence of the surface.

    Arguments
    ---------
    model : 'line' or 'sine'
        Model the discplacement sensitivity as line or as a modulated sine.
    fun_kw : None or dict with keyword arguments for the oscillation function.
    """
    h = height
    r = radius

    h0 = parameters['h0'].value
    fs = parameters['focal_shift'].value
    #h0_ = parameters['h0_'].value

    lin = _lin_fun(h, x0=h0, y0=dissens_surf, slope=fs * slope)
    rel_drag = _rel_drag_fit_fun(parameters, h, r,
                                 lateral=lateral, method=method)

    if fit_osci:
        A = fun_kw['A']
        d = fun_kw['d']
        phi = fun_kw['phi']
        n = fun_kw['ref_ind']
        wl = fun_kw['wavelength']
        osci = A * exp(-d * h) * sin(4 * pi / (wl / n) * fs * h + phi)
        #osci = (A + d * h) * sin(4 * pi / (wl / n) * fs * h + phi)
        #osci = A * exp(-(d * (h+h0_))**2 * sin(4 * pi / (wl / n) * fs * h + phi))

        return lin * rel_drag**0.5 + osci
    else:
        return lin * rel_drag**0.5


def _kappa_fit_fun(parameters,
                   height,  # in um
                   trap_stiffness_surf,  # in pN/nm
                   slope,  # in pN/(nm*um)
                   radius,  # in um
                   lateral=True,
                   method='viscosity'
                   ):
    """
    Return the value(s) of the theoretically acquired experimental
    trap_stiffness, at a given height above the surface and
    neglecting the presence of the surface.
    """
    h = height
    r = radius

    h0 = parameters['h0'].value
    fs = parameters['focal_shift'].value

    lin = _lin_fun(h, x0=h0, y0=trap_stiffness_surf, slope=fs * slope)
    rel_drag = _rel_drag_fit_fun(parameters, h, r,
                                 lateral=lateral, method=method)
    return lin / rel_drag


def _height_fit_fun(parameters,
                    names,
                    heights_dict,  # in um
                    radius,  # in um
                    data_dict=None,  # in units [1, nm/mV, pN/nm]
                    errors_dict=None,  # in units [1, nm/mV, pN/nm]
                    lateral_dict=None,
                    method='viscosity',
                    return_fit_dict=False,
                    fit_dissens_osci=False,
                    dissens_fun_kw={}):
    """
    Calculate the residuals of the height-dependent fit.

    If no data is given the function returns the evaluated function for the
    given parameters.

    Arguments
    ---------
    parameters : lmfit.Parameters
        Parameters object holding the lmfit.Parameter objects: beta_, mbeta_,
        kappa_, mkappa_, focal_shift, corr and h0, where "_" referes to the
        names of the axes.
    names : list of str
        List holding the names of the axes.
    height_dict : OrderedDict
        Dictionary holding the height vectors for the different axes and the
        height vector for the "rel_drag" to be fitted. So the following key
        names are expected: 'kappa_x_pc', 'beta_x_pc' etc.
    radius : float
        Specified radius of the microsphere in micrometers
    focal_shift : float
    data_dict : OrderedDict
        Dictionary of the data to be fitted with same keys as height_dict.
    errors_dict : OrderedDict
        Dictionary of the errors of the data to be fitted with same keys as
        height_dict.
    lateral : OrderedDict
        Dictionary with same keys as height_dict refering to whether the
        axis is in the lateral or axial direction. E.g. for 'x', 'y', 'z' and
        'rel_drag' as names for the three axes and the drag:
         lateral_dict = {'x': True, 'y': True, 'z': False}
    method : str
        Either 'viscosity' or 'radius', referes to the fitting routine used,
        where the drag is either calculated by:
         drag = 6 * pi * corr * eta * r * Faxen(r)
        or
         drag = 6 * pi * eta * corr * r * Faxen(corr * r).
    return_fit_dict : bool
        Return a dictionary (True) or a flattened array (False)
    fit_dissens_osci : bool
        Whether to try to fit oscillations on the displacement sensitivity. If
        True, the fit will have an additional summand of the shape:
        A * exp(d * heights) * sin( 4 * pi / (n * wavelength) * focal_shift *
        heigths + phi).
    dissens_fun_kw : dict
        If fit_dissens_osci is True, this dictionary needs to be provided. It
        holds the following key-value pairs:

         - ref_ind -- refractive index
         - wavelength -- wavelegth of the detection (trapping) laser in vacuum


    Returns
    -------
    residuals : array or dict
    """
    p = parameters.valuesdict()

    beta_ = OrderedDict()    # value at true surface
    mbeta_ = OrderedDict()   # slope with respect to apparent distance
    kappa_ = OrderedDict()   # value at true surface
    mkappa_ = OrderedDict()  # slope with respect to apparent distance

    for name in names:
        beta_[name] = p['beta_' + name]
        mbeta_[name] = p['mbeta_' + name]
        kappa_[name] = p['kappa_' + name]
        mkappa_[name] = p['mkappa_' + name]

    if fit_dissens_osci:
        A_ = OrderedDict()
        d_ = OrderedDict()
        phi_ = OrderedDict()
        for name in names:
            A_[name] = p['A_' + name]
            d_[name] = p['d_' + name]
            phi_[name] = p['phi_' + name]

    rel_drag = _rel_drag_fit_fun(parameters,
                                 heights_dict['rel_drag'],
                                 radius,
                                 lateral=lateral_dict['rel_drag'],
                                 method=method)

    model_ = [list(rel_drag)]
    model_dict = OrderedDict()
    model_dict['rel_drag'] = rel_drag

    for name in names:
        if fit_dissens_osci:
            dissens_fun_kw['A'] = A_[name]
            dissens_fun_kw['d'] = d_[name]
            dissens_fun_kw['phi'] = phi_[name]

        dissens = _dissens_fit_fun(parameters,
                                   heights_dict['beta_' + name + '_pc'],
                                   beta_[name],
                                   mbeta_[name],
                                   radius,
                                   lateral=lateral_dict[name],
                                   method=method,
                                   fit_osci=fit_dissens_osci,
                                   fun_kw=dissens_fun_kw)
        model_.append(dissens)
        model_dict['beta_' + name + '_pc'] = dissens

        kappa = _kappa_fit_fun(parameters,
                               heights_dict['kappa_' + name + '_pc'],
                               kappa_[name],
                               mkappa_[name],
                               radius,
                               lateral=lateral_dict[name],
                               method=method)
        model_.append(kappa)
        model_dict['kappa_' + name + '_pc'] = kappa

    result_flat = array(list(flatten_list(model_)))

    if not return_fit_dict:
        if data_dict is None:
            return result_flat

        data_ = [data_dict['rel_drag']]
        for name in names:
            data_.append(data_dict['beta_' + name + '_pc'])
            data_.append(data_dict['kappa_' + name + '_pc'])
        data_flat = array(list(flatten_list(data_)))

        if errors_dict is None:
            return (data_flat - result_flat)
        else:
            err_ = [errors_dict['rel_drag']]
            for name in names:
                err_.append(errors_dict['beta_' + name + '_pc'])
                err_.append(errors_dict['kappa_' + name + '_pc'])
            err_flat = array(list(flatten_list(err_)))
            return ((data_flat - result_flat) / err_flat)
    else:
        if data_dict is None:
            return model_dict

        result = OrderedDict()
        if errors_dict is None:
            for k, data in data_dict.items():
                result[k] = data - model_dict[k]
            return result
        else:
            for k, data in data_dict.items():
                result[k] = ((data - model_dict[k]) / errors_dict[k])
            return result


def fit_height_data(names,
                    heights_dict,  # in um
                    data_dict,  # in units [1, nm/mV, pN/nm]
                    errors_dict,
                    radius,  # in um
                    parameters=None,
                    fit_drag_first=True,
                    drag_fit_method='radius',  # viscosity or radius
                    lateral_dict=None,
                    fit_dissens_osci=False,
                    fit_fun_kw={},
                    fit_focal_shift=False,
                    iterations=1,
                    fit_report=False,
                    verbose=False,
                    **init_kwargs):
    """
    Fit the height dependent values for drag, displacement sensitivity and
    trap_stiffness globally.

    Arguments
    ---------
    names : list of str
        List holding the names of the axes.
    height_dict : OrderedDict
        Dictionary holding the height vectors for the different axes and the
        height vector for the "rel_drag" to be fitted. So the following key
        names are expected: 'kappa_x_pc', 'beta_x_pc' etc.
    data_dict : OrderedDict
        Dictionary of the data to be fitted with same keys as height_dict.
    errors_dict : OrderedDict
        Dictionary of the errors of the data to be fitted with same keys as
        height_dict.
    radius : float
        Specified radius of the microsphere in micrometers.
    focal_shift : float
        Factor that describes the physical shift of the trap above the surface.
    parameters : lmfit.Parameters
        Parameters object holding the lmfit.Parameter objects: beta_, mbeta_,
        kappa_, mkappa_, corr and h0, where "_" referes to the names of the
        axes. See gen_height_fit_pars().
    fit_drag_first : bool
        Whether to fit the relative drag independently first.
    drag_fit_method : {'radius', 'viscosity'}
        Determines if the correction factor 'corr' is associated with the
        'radius' or the 'viscosity'.
    lateral_dict : dict
        dictionary that hold the information, if the axis with a certain name
        (the key) is considered to be lateral (True), i.e. the movement is
        parallel to the surface, or axial (False).
    fit_dissens_osci : bool
        Whether to try to fit oscillations on the displacement sensitivity. If
        True, the fit will have an additional summand of the shape:
        A * exp(d * heights) * sin( 4 * pi / (n * wavelength) * focal_shift *
        heigths + phi).
    fit_fun_kw : dict
        If fit_dissens_osci is True, this dictionary needs to be provided. It
        holds the following key-value pairs:

         - ref_ind -- refractive index
         - wavelength -- wavelegth of the detection (trapping) laser in vacuum

    iterations : int
        Number of iteration. This might not have an affect to the result.
    verbose : bool
        Give more info.

    Keyword Arguments
    -----------------
    kwargs
        Provide initial guess values for the fit by setting the value of the
        parameter. See gen_height_fit_pars() for possible values.

    Returns
    -------
    results : tuple
        A 2-tuple of the minimizer results of the global fit and the fit to
        the relative drag if *fit_rel_drag* is True.

    Note
    ----
    The idea behind the drag_fit_method is to account for possible deviations
    from the 'known' viscosity or the radius, which would result in a false
    normalization of the drag to the stokes drag.
    I.e.:
        drag = drag_stokes' * faxen's law
             = 6 * pi * viscosity' * radius' * FL(radius', height)
             = 6 * pi * c_visc * viscosity * c_rad * radius
                 * FL(c_rad * radius, height)
             = drag_stokes * c_visc * c_rad * FL(c_rad * radius, height)
    where viscosity (of water) can be calculated if the temperature is known
    and radius is known from the bead manufacturer. Thus, one can write:
        rel_drag = drag / drag_stokes
                 = c_visc * c_rad * FL(c_rad*radius, height)

    The idea is now to set one of the two correction factors fixed to 1.0 and
    let the other factor being varied within the fitting algorithm. The method
    specifies which correction factor is allowed to vary.
    """
    if parameters is None:
        if 'corr' not in init_kwargs:
            init_kwargs['corr'] = min(data_dict['rel_drag'])

        pn = {}
        if fit_dissens_osci:
            pn = {'A': (1, -10, 10),
                  'd': (0, -10, 10),
                  'phi': (0, -pi/2, pi/2)}

        pars = gen_height_fit_pars(names=names,
                                   pars_named=pn,
                                   **init_kwargs)
    else:
        pars = parameters

    if fit_focal_shift:
        pars['focal_shift'].vary = True
    else:
        pars['focal_shift'].vary = False

    if verbose:
        print('Drag fit method: {}'.format(drag_fit_method))

    if lateral_dict is None:
        lateral_dict = {}
        lateral_dict['rel_drag'] = True
        for name in names:
            lateral_dict[name] = True
        if verbose:
            print('All axes are condsidered parallel to the surface.')

    if fit_drag_first:
        drag_minimizer = fit_rel_drag(heights_dict['rel_drag'],
                                      data_dict['rel_drag'],
                                      errors_dict['rel_drag'],
                                      radius,
                                      parameters=pars,
                                      lateral=lateral_dict['rel_drag'],
                                      method=drag_fit_method,
                                      verbose=verbose)
    else:
        drag_minimizer = None

    if verbose:
        print('Global fit(s):')

    for i in range(iterations):
        corr = pars['corr'].value
        height_minimizer = minimize(_height_fit_fun,
                                    pars,
                                    args=(names,
                                          heights_dict,
                                          radius),
                                    kws={'data_dict': data_dict,
                                         'errors_dict': errors_dict,
                                         'lateral_dict': lateral_dict,
                                         'method': drag_fit_method,
                                         'return_fit_dict': False,
                                         'fit_dissens_osci': fit_dissens_osci,
                                         'dissens_fun_kw': fit_fun_kw},
                                    nan_policy='propagate')
        if verbose:
            print('Iteration step {0}: Correction factor '
                  ' = {1:1.5f}'.format(i, corr))
            if i == iterations-1:
                report_fit(height_minimizer, show_correl=True)
            else:
                report_fit(height_minimizer, show_correl=False)

        for k, param in pars.items():
            if (param.value == param.min):
                if verbose:
                    print('Parameter "{0}" hit the lower bound at {1:1.3f}.'
                          'The bound was decreased by 30% to {2:1.3f}'
                          ''.format(k, param.min, 0.7*param.min))
                param.min *= 0.7
            elif (param.value == param.max):
                if verbose:
                    print('Parameter "{0}" hit the upper bound at {1:1.3f}.'
                          'The bound was increased by 30% to {2:1.3f}'
                          ''.format(k, param.min, 0.7*param.max))
                param.max *= 1.7

    if height_minimizer.success:
        if verbose:
            print('Least squares fit to the height dependent data succeeded.')
            report_fit(height_minimizer)

    else:
        if verbose:
            print('Least squares fit to the height dependent data did not '
                  'succeed.')
            print('Minimizer message: "{0:s}"'
                  ''.format(height_minimizer.message))
            print('lmdif message: "{0:s}"'
                  ''.format(height_minimizer.lmdif_message))
    # print fit results anyway
    if fit_report and verbose is False:
        report_fit(height_minimizer, show_correl=False)

    return (height_minimizer, drag_minimizer)


class HeightFitResult(object):
    """
    Height fit result.

    You can access the height fit results by self.'paramname' or
    self.'paramname_err' to get, e.g. 'beta_x' or 'mkappa_z_err'.
    """
    def __init__(self,
                 method='',
                 minimizer=None):
        self.method = method
        self.minimizer = minimizer

    def eval(self, name):
        return self.minimizer.eval[name]

    def residuals(self, name):
        return self.minimizer.residuals[name]

    def get_params(self):
        return self.minimizer.params

    @property
    def params(self):
        return self.get_params()

    def get_h0(self, unit='um', error=False):
        if unit == 'um':
            conv = 1.0
        else:
            conv = ureg('um').to(unit).magnitude

        if error:
            output = self.params['h0'].stderr * conv
        else:
            output = self.params['h0'].value * conv

        return output

    @property
    def h0(self):
        return self.get_h0()

    @property
    def h0_err(self):
        return self.get_h0(error=True)

    def get_dissens(self, name, unit='nm/mV', error=False, slope=False):
        if unit == 'nm/mV':
            conv = 1.0
        else:
            conv = ureg('nm/mV').to(unit).magnitude

        if slope:
            key = 'm'
        else:
            key = ''

        if error:
            return self.params[key + 'beta_' + name].stderr * conv
        else:
            return self.params[key + 'beta_' + name].value * conv

    def get_trap_stiffness(self, name, unit='pN/nm', error=False, slope=False):
        if unit == 'pN/nm':
            conv = 1.0
        else:
            conv = ureg('pN/nm').to(unit).magnitude

        if slope:
            key = 'm'
        else:
            key = ''

        if error:
            return self.params[key + 'kappa_' + name].stderr * conv
        else:
            return self.params[key + 'kappa_' + name].value * conv

    def __getattr__(self, name):
        """
        """
        if name.startswith('m'):
            slope = True
            name_ = name[1:]
        else:
            slope = False
            name_ = name

        if name_.startswith('beta'):
            if name_.endswith('err'):
                error = True
            else:
                error = False

            axis = name_.split('_')[1]
            return self.get_dissens(axis, error=error, slope=slope)

        elif name.startswith('kappa'):
            if name_.endswith('err'):
                error = True
            else:
                error = False

            axisname = name_.split('_')[1]
            return self.get_trap_stiffness(axisname, error=error, slope=slope)
        else:
            raise AttributeError(name)


class HeightCalibration(object):
    """
    Manage height-dependent psd measurements and their fitting.
    """
    def __init__(self):
        self._mask = None
        self._outlier_mask = {}
        self._names = []
        self._laterality = {}
        self._heights = []
        self._height_unit = 'um'
        self.psdfits = OrderedDict()
        self._arrays = OrderedDict()

        self._dissens = OrderedDict()
        self._dissens_unit = 'nm/mV'
        self._trap_stiffness = OrderedDict()
        self._trap_stiffness_unit = 'pN/nm'
        self._drag = OrderedDict()
        self._drag_unit = 'nN*s/m'
        self._red_chi2 = OrderedDict()
        self._offset = OrderedDict()
        self._temp = array([])
        self._temp_unit = 'K'

        self._radius = None
        self._radius_err = None
        self._radius_unit = 'um'
        self._ex_axes = None

        self._wl = 1064
        self._wl_unit = 'nm'
        self.focal_shift = 0.800
        self.focal_shift_err = 0.001
        self.ref_ind = 1.326  # refractive index of water @ 25 °C
        # maximum relative drag used for fitting routine
        self.max_rel_drag = 1.7
        self.height_range = (-inf, inf)
        self.conf_level = 0.95  # used for outlier determination
        self.fit_psd_kwargs = {}
        self.height_fit = None
        self.height_fit_results = None
        self.rel_drag_fit = None

        self.height_offset = 0.0

        self.directory = None
        self.basename = None

        self.datafilename = None
        self.paramfilename = None

    def __getattr__(self, name):
        """
        Shortcut for HeightCalibration().arrays['name'].
        """
        if name in self._arrays.keys():
            return self._arrays[name]
        else:
            raise AttributeError(name)

    @property
    def names(self):
        return self._names

    def get_heights(self,
                    name=None,
                    unit=None,
                    getQuantity=False,
                    nomask=False,
                    get_outliers=False,
                    inverse_mask=False
                    ):
        """
        Return the array of heights in units specified by unit.

        If unit is None self._height_unit is used.

        Arguments
        ---------
        unit : str
            Unit of dimension length.
        name : str
            Name of the axis. passed over to self.get_mask().
        getQuantity : bool
            Whether to return a pint.Quantity (True) or an numpy array.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name=name,
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)

        if unit is None or unit == self._height_unit:
            unit = self._height_unit
            output = array(self._heights)[~mask]
        else:
            conv = ureg(self._height_unit).to(unit).magnitude
            output = array(self._heights)[~mask] * conv
        if getQuantity:
            return output * ureg(unit)
        else:
            return output

    def add_height_offset(self, offset, unit='um'):
        """
        Add an offset to all heights.

        This primarily, is useful if psdfits are present and the heights need
        to be adjusted, to ensure 'correct' height is used for fitting.
        """
        offset_ = offset * ureg(unit).to(self._height_unit).magnitude
        heights = self._heights  # list of heights

        h_ = []

        if self.psdfits:
            pf_ = OrderedDict()

            for height in heights:
                h = height + offset_
                h_.append(h)

                pf = self.psdfits[height]
                pf.psdm.exp_setting.set_height(h, self._height_unit)
                pf_[h] = pf
            self.psdfits = pf_
        else:
            print('Only heights were adjusted. PSDFit objects are not'
                  ' present. This is probably due to loading of the height'
                  ' dependent data from a file.')
        self._heights = h_ or [height + offset_ for height in heights]
        self.height_offset += offset_

    def set_height_unit(self, unit, convert=True):
        """
        Set the unit of the height vector. If convert is True the heights are
        converted into this unit, if False their magnitude is not changed, but
        the unit.
        """
        unit_ = str(ureg(unit).units)
        if convert:
            conv = ureg(self._height_unit).to(unit).magnitude
            self._heights *= conv
        self._height_unit = unit_

    @property
    def radius(self):
        if self._radius is not None:
            return (self._radius * ureg(self._radius_unit))
        else:
            return None

    def get_radius(self, unit=None):
        """
        Return the radius in units specified by unit.

        If unit is None radius is returned in units given by self._radius_unit.
        """
        if unit is None:
            return self._radius
        else:
            conv = ureg(self._radius_unit).to(unit).magnitude
            return self._radius * conv

    def set_radius(self, radius, unit='um'):
        """
        Set the radius of the used bead in units 'unit'.
        """
        self._radius = radius
        if self._radius_err is not None:
            if self._radius_unit is not str(ureg(unit).units):
                conv = ureg(self._radius_unit).to(unit).magnitude
                self._radius_err *= conv

        self._radius_unit = str(ureg(unit).units)

    @property
    def radius_err(self):
        if self._radius_err is not None:
            return self._radius_err * ureg(self._radius_unit)
        else:
            return None

    def get_radius_err(self, unit=None):
        """
        Return the radius in units specified by unit.
        """
        if unit is None:
            return self._radius_err
        else:
            conv = ureg(self._radius_unit).to(unit).magnitude
            return self._radius_err * conv

    def set_radius_err(self, error, unit='um'):
        """
        Set the error of the radius.

        If the unit is not the same as the unit of the radius, the error is
        converted into the same unit.
        """
        if self._radius_unit is not str(ureg(unit).units):
            conv = ureg(unit).to(self._radius_unit).magnitude
            error *= conv
        self._radius_err = error

    @property
    def wavelength(self):
        return self.get_wavelength()

    def get_wavelength(self, unit=None):
        """ Return the wavelength in units specified by unit. """
        if unit is None:
            return self._wl
        else:
            conv = ureg(self._wl_unit).to(unit).magnitude
            return self._wl * conv

    def set_wavelength(self, wavelength, unit='nm'):
        self._wl = wavelength
        self._wl_unit = str(ureg(unit).units)

    def reset_mask(self, outliers=False, conf_level=None):
        """
        Reset the mask. If outliers is True also reset the outliers mask.

        Arguments
        ---------
        outliers : bool
            Whether to reset the outliers masks, as well. This option runs
            generate_ouliers_mask() method.
        conf_level : None or float between 0 and 1
            If None, self.conf_level is used as confidence level for the
            outliers mask. If a value between 0 and 1 is given, self.conf_level
            is set to that value and the outliers are generated with that value.
        """
        self._mask = zeros(shape(self._heights)) > 0
        if outliers:
            self.generate_outliers_masks(conf_level=conf_level)

    def get_mask(self, name=None, get_outliers=False, inverse=False):
        """
        Return mask of the given axis. If axis is None the outlier mask is
        omitted.

        Arguments
        ---------
        name : str
            Name of the axis, e.g. 'x' or 'ex_axis' which would account for
            changing excitation axes with height
        """
        if self._mask is None:
            self.reset_mask()
        if name is None:
            return self._mask
        elif name == 'ex_axis':
            names = self._names
            ex_axes = self._ex_axes
            mask_ = []
            for idx, axis_number in enumerate(ex_axes):
                name_ = names[int(axis_number)]
                mask_.append(self._outlier_mask[name_][idx])
            mask = array(mask_)
            if get_outliers:
                mask = ~mask
            output = logical_or(self._mask, mask)
        else:
            mask = self._outlier_mask[name]
            if get_outliers:
                mask = ~mask
            output = logical_or(self._mask, mask)

        if inverse:
            return ~output
        else:
            return output

    def exclude_height(self, height, unit='um', reset=False):
        """
        Exclude values for height 'height'.

        Arguments
        ---------
        height : float
        unit : str
            Unit of height.
        reset : bool
            If True the internal mask is reset (the outlier mask is not
            affected) - see reset_mask() for this purpose.
        """
        if reset:
            self.reset_mask(outliers=False)

        heights = array(self._heights) * ureg(self._height_unit)
        if not(isinstance(height, ureg.Quantity)):
            height = height * ureg(unit)
        mask = height == heights
        self._mask = logical_or(self._mask, mask)

    def exclude_by_max_rel_drag(self, max_rel_drag=None, reset=False):
        """
        Exclude values for relative drags above max_rel_drag.

        If max_rel_drag is None the internal self.max_rel_drag is used.

        Arguments
        ---------
        max_rel_drag : float or Nona
        reset : bool
            If True the internal mask is reset (the outlier mask is not
            affected) - see reset_mask() for this purpose.
        """
        if reset:
            self.reset_mask(outliers=False)
        if max_rel_drag is None:
            max_rel_drag = self.max_rel_drag
        else:
            self.max_rel_drag = max_rel_drag

        rel_drag = self.get_drag(relative_drag=True, nomask=True)

        mask = rel_drag > max_rel_drag
        self._mask = logical_or(self._mask, mask)

    def exclude_heights_outside(self, hmin, hmax, unit='um', reset=False):
        """
        Exclude values for heights outside [hmin, hmax].

        Arguments
        ---------
        hmin : float
        hmax : float
        unit : str
            Unit of hmin and hmax.
        reset : bool
            If True the internal mask is reset (the outlier mask is not
            affected) - see reset_mask() for this purpose.
        """
        if reset:
            self.reset_mask(outliers=False)

        if unit != self._height_unit:
            hmin *= ureg(unit).to(self._height_unit).magnitude
            hmax *= ureg(unit).to(self._height_unit).magnitude

        heights = self.get_heights(unit=unit, nomask=True)
        mask = logical_or(heights < hmin, heights > hmax)
        self._mask = logical_or(self._mask, mask)

    def generate_outliers_masks(self, names=None, conf_level=None):
        """
        Generate masks for all available axes to mask outliers.

        Outliers are - for a given level of confidence - identified by the
        reduced chi² value and the number of degrees of freedom of the
        fit to the power spectral density. 

        Arguments
        ---------
        names : str or list of str
            Axis name(s) to generate the outliers mask(s) for.
        conf_level : None or float between 0 and 1
            If None, self.conf_level is used as confidence level for the
            outliers mask. If a value between 0 and 1 is given,
            the outliers are generated according to that.

        See also
        --------
        is_outlier()
        """
        if conf_level is None:
            conf_level = self.conf_level

        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self._names:
            redch2 = self._red_chi2['red_chi2_' + name]
            nfree = self._red_chi2['N_free_' + name]
            mask = array(is_outlier(redch2, nfree, conf_level))
            self._outlier_mask[name] = mask

    def reset_recalc(self):
        """
        Reset all calculated data.

        This can be useful after, for example, adjusting the height by an
        offset. The method clears, all data taken from psdfits and initiates a
        re-collection of the data from the psd fits if. If psdfits are not
        present the height-mask is reset, only.
        """
        self.reset_mask(outliers=True)

        if self.psdfits:
            self._arrays = OrderedDict()
            self._dissens = OrderedDict()
            self._trap_stiffness = OrderedDict()
            self._drag = OrderedDict()
            self._red_chi2 = OrderedDict()
            self._offset = OrderedDict()
            self._gen_height_data()
            self._construct_arrays()
        

    def get_values(self, vname, name=None, nomask=False):
        """
        Return the values specified by name.

        If axis is not specified the outlier mask is omitted.

        Arguments
        ---------
        vname : str
            The name of the values that shall be returned.
        name : str
            Specify the name of the axis, e.g. 'x'. The returned vector is
            then masked with the outlier mask, as well.
        """
        if len(self._arrays) <= 0:
            self._construct_arrays()
        if nomask:
            return self._arrays[vname]
        else:
            mask = self.get_mask(name=name)
            return self._arrays[vname][~mask]

    def reset_psd_masks(self):
        """ Reset the psd masks. """
        for psdfit in self.psdfits.values():
            psdfit.psdm.reset_masks()

    #####----------------------------------------------------------------------
    #---- manage psd fits  ----
    #####----------------------------------------------------------------------

    def add_psdfit(self, height, psdfit, unit='um', round_to=3):
        """
        Add a **PSDFit** object for a **height** in the specified unit.

        The psdfit objects are stored in the psdfits dictionary.

        Arguments
        ---------
        height : float
            Height of the psd measurement.
        psdfit : PSDFit
            PSDFit object.
        unit : str
        round_to : int
            Digit to which the height value is rounded to.
        """
        if not(isinstance(psdfit, PSDFit)):
            raise Exception('You should add PSDFit objects.')

        height_ = (height * ureg(unit)).to(self._height_unit).magnitude
        height = round(height_, round_to)

        self._heights.append(height)
        self.psdfits[height] = psdfit

        r = psdfit.psdm.exp_setting.get_radius(unit=self._radius_unit)
        if self.radius is None:
            self.set_radius(r, self._radius_unit)
        else:
            ## TODO double check units before comaring magnitudes
            if len(set([self.get_radius(), r])) > 1:
                raise Exception('Bead size must be equal. Given radii '
                                'are: {0:1.3~} and {1:1.3~}'
                                ''.format(self.radius, r))
        dr = psdfit.psdm.exp_setting.get_radius_err(unit=self._radius_unit)
        self.set_radius_err(dr, unit=self._radius_unit)

        if self._names == []:
            self._names = psdfit.names
        else:
            if len(set([''.join(self._names), ''.join(psdfit.names)])) > 1:
                raise Exception('Axis names are not equal: got {0} and {1} at '
                                ' height {2:1.3e} {3:s}'
                                ''.format(self._names, psdfit.names, height,
                                          unit))
        if self._laterality == {}:
            self._laterality = psdfit.psdm.get_laterality()

    def sort_by_height(self):
        """
        Sort the *PSDFit* objects with respect to acending height.
        """
        d = OrderedDict()
        self._heights.sort()
        for height in self._heights:
            d[height] = self.psdfits[height]
        self.psdfits = d

    def get_psd_files(self,
                      basename,
                      directory,
                      datafile_ext='.dat',
                      paramfile_suffix='_psd_parameters',
                      paramfile_ext='.txt',
                      verbose=False):
        """
        Scan the given directory for psd-data- and parameter-files.

        The files are collected in self.files.

        Arguments
        ---------
        basename : str
            String that all height calibration files have in common at the
            beginning of the file name.
        directory : path
            Path to the directory.
        datafile_ext : str
            Extension of the data file.
        paramfile_suffix : str
            Suffix that defines the parameter file.
        paramfile_ext : str
            Parameter file extension.

        Note
        ----
        The files can be accessed through the 'files' Attribute.

        """
        self.basename = basename
        self.directory = directory

        datafiles = [f for f in listdir(directory)
                     if (f.startswith(basename) and f.endswith(datafile_ext))]

        paramfile_end = ''.join([paramfile_suffix, paramfile_ext])

        paramfiles = [f for f in listdir(directory)
                      if (f.startswith(basename) and
                          f.endswith(paramfile_end))]

        files = []
        for df in datafiles:
            pf = ''.join([df.split(datafile_ext)[0], paramfile_end])
            if pf not in paramfiles:
                print('"{}" omitted - parameter file not found.'.format(df))
            else:
                if verbose:
                    print('found: {}'.format(df))
                files.append((df, pf))

        files.sort()

        self.files = files

    def init_psdfit(self,
                    datafile,
                    paramfile,
                    directory=None,
                    bounds=None,
                    f_exclude=None,
                    verbose=False,
                    **kwargs):
        """
        Generate a psdfit object of the psd measurement,for the given
        data- and parameter file.

        Arguments
        ---------
        datafile : complete name of the data file.
        paramfile : complete name of the parameter file.
        directory : path to the data directory. If None, the attribute
            directory is checked.
        bounds : tuple (fmin, fmax) or dict {name: (fmin, fmax), ...}
            Set the bound for all (list) or specific names (dict).
        f_exclude : list or dict
            Specify which datapoints to exclude. If a list is given, the
            frequencies are excluded for all names. If a dict is given the
            keys should correspond to the names of the axes and their values
            can be either lists or floats.

        Keyword Arguments
        -----------------
        keyword-arguments passed to *PSDFit* and their to the
        analytical_lorentzian_fit method.

        See also:
        ---------
        PSDFit
        """
        if directory is None:
            if self.directory:
                directory = self.directory
            else:
                raise Exception('No directory given.')

        psdm = PSDMeasurement(warn=False)
        psdm.load(directory, datafile, paramfile=paramfile)

        height = psdm.exp_setting.get_height(unit='nm')

        if verbose:
            print('"{0:s}" loaded for height {1:1.1f} nm'
                  ''.format(datafile, height))

        # exclude specified frequencies
        if f_exclude:
            # check if axis names are specified or if frequencies for all
            # axes schould be excluded
            if hasattr(f_exclude, 'keys'):
                for name, freqs in f_exclude.items():
                    psdm.exclude_freq(freqs, names=name)
            else:
                psdm.exclude_freq(f_exclude)

        pf = PSDFit(psdm, bounds=bounds, verbose=verbose, **kwargs)
        self.add_psdfit(height, pf, unit='nm')

    def gen_psd_fits(self,
                     bounds=None,
                     f_exclude=None,
                     verbose=False,
                     **kwargs):
        """
        Initialize a psdfit object for each data and parameter file in the
        'files' attribute.

        Boundaries and frequencies that should be excluded can be provided,
        thus they are also not considered during the initialization of the
        psdfit object.

        Arguments
        ---------
        bounds : tuple (fmin, fmax) or dict {name: (fmin, fmax), ...}
            Set the bound for all (list) or specific names (dict).
        f_exclude : list or dict
            Specify which datapoints to exclude. If a list is given, the
            frequencies are excluded for all names. If a dict is given the
            keys should correspond to the names of the axes and their values
            can be either lists or floats.
        verbose : bool
            Give more info.

        Keyword Arguments
        -----------------
        keyword-arguments passed to *PSDFit* and their to the
        analytical_lorentzian_fit method.

        See also:
        ---------
        get_psd_files, init_psdfit

        """
        if not self.files:
            raise Exception('No files found. Run get_psd_files() first.'
                            'Alternatively, you can also set the "files" '
                            'Attribute with a list of (datafile, paramfile) '
                            'tuples.')

        if verbose:
            print('generating PSDFit objects and '
                  'adding them to self.psdfits')
        else:
            print('Loading data and preparing fits to PSDs.')

        for dfile, pfile in self.files:
            if verbose:
                print('loading {} with {}'.format(dfile, pfile))

            self.init_psdfit(dfile,
                             pfile,
                             bounds=bounds,
                             f_exclude=f_exclude,
                             verbose=verbose,
                             **kwargs)
            if not verbose:
                print('.', flush=True, end='')

        self.reset_mask()
        print('')

    #####----------------------------------------------------------------------
    #---- run psd fits ----
    #####----------------------------------------------------------------------

    def setup_psd_fits(self,
                       height=None,
                       names=None,
                       model='hydro',
                       lp_filter=False,
                       lp_fixed=False,
                       aliasing=False,
                       f_sample=None,
                       N_alias=9,
                       debias=True,
                       f3dB=None,
                       alpha=None,
                       dynamic=False,
                       dyn_factor=1.0,
                       **fit_psd_kws):
        """
        Setup the psd fits, i.e. the model to use and whether to consider a
        low-pass filter or aliasing of the data.

        Arguments
        ---------
        height : float or None
            height of the fit that should be set up.
        model : {'hydro', 'lorentzian'}
            Whether to use the hydrodynamically correct or a lorentzian-shaped
            psd as a model function.
        lp_filter : bool
            Add a low-pass filter to the model function.
        lp_fixed : bool
            Whether the low-pass filter parameters are fixed (True) or can be
            varied during the fit (False).
        aliasing : bool
            Consider aliasing in the model function. If True, f_sample needs to
            be provided.
        f_sample : float or None
            If aliasing is True, the sampling frequency of the data must be
            given.
        N_alias : int
            At which iteration to stop the infinite sum (see aliasing). The
            default (N_alias=9) give reasonable results (<0.5 % deviation).
        debias : bool
            Whether to account for biasing of the least-squares fitting.
            This changes the diffusion constant by a factor of n/(n+1).
        f3dB : float
            The cut-off frequency of the low-pass filter.
        alpha : float
            The amount of the signal which is not low-pass filtered.
        dynamic : bool
            Whether to determine the model to be used dynamically, depending
            on the given height. Below the threshold a lorentzian is used.
        dyn_factor : float
            Factor, multiplied by the radius of the trapped particle, that
            determines the distance at which to switch the model function.
        **fit_psd_kws : keyword arguments passed to PSDFit.setup_fit().

        Notes
        -----
        When using the dynamic option, consider using add_height_offset()
        first.
        """
        thrsh = dyn_factor * self.get_radius(unit=self._height_unit)

        if height:
            try:
                heights = iter(height)
            except:
                heights = [height]
        else:
            heights = self._heights

        for height in heights:
            if dynamic:
                model = 'hydro' if height > thrsh else 'lorentzian'

            self.psdfits[height].setup_fit(names=names,
                                           model=model,
                                           lp_filter=lp_filter,
                                           lp_fixed=lp_fixed,
                                           aliasing=aliasing,
                                           f_sample=f_sample,
                                           N_alias=N_alias,
                                           debias=debias,
                                           f3dB=f3dB,
                                           alpha=alpha,
                                           **fit_psd_kws)

    def fit_psds(self,
                 heights=None,
                 names=None,
                 dynamic_bounds=False,
                 bounds=None,
                 f_exclude=None,
                 use_heights=False,
                 fitreport=False,
                 plot_fits=False,
                 plot_delay=0.0,
                 plot_axis=None,
                 figsize=(9, 6),
                 save_plots=False,
                 save_as='png',
                 directory=None,
                 basename=None,
                 verbose=False,
                 plot_kws={},
                 **kwargs):
        """
        Fit the psds of the psdfit object at the corresponding heights.

        If No height is provided all are fitted.

        Arguments
        ---------
        heights : float or list of floats
            Height or Heights to run the psd fit of.
        names : str or list of strings
            Names of the axes to be fitted.
        dynamic_bounds : bool
            Whether to determine the fit boundaries from the corner frequency
            that the analytical fit resulted in. This would lead to
            dynamically changing boundaries.
        bounds : tuple (fmin, fmax) or dict {'name': (fmin, fmax), ...}
            defines the lower and upper frequency bound. All data outside the
            interval is masked.
        f_exclude : float or list of floats or dictionary that specifies the
            frequencies for each name. 'names' then has no effect.
        use_heights : bool
            Whether to consider the height when using the hydrodynamically
            correct psd as a model function. Infinity is used if set to False.
        fitreport : bool
            Whether to rpint a fit report.
        plot_fits : bool
            Plot the fits.
        plot_delay : float
            Whether to wait for this amount of seconds until the next fit is
            done.
        plot_axis : matplotlib.Axis
            Axis to plot the results to.
        figsize : tuple(width, height)
            Size of the figure.
        save_plots : bool
            Save the plot.
        save_as : str
            Format of the save image, e.g. 'png', 'svg', etc.
        directory : str
            Where to store the images
        basename : str
            Basename of the images.
        verbose : bool
            Give more info.

        Keyword Arguments
        -----------------
        Passed over to PSDFit.fit_psd(), e.g. initial guess values for 'D' and
        'f_c'.
        """
        if save_plots:
            conv = ureg(self._height_unit).to('nm').magnitude
            if directory is None:
                directory = self.directory if self.directory else './'
            if basename is None:
                bn = self.basename if self.basename else 'psd_fits'
        if plot_fits:
            plt.ion()
            fig = plt.figure(figsize=figsize)
            if not plot_axis:
                plot_axis = fig.gca()

        if heights:
            try:
                heights_ = iter(heights)
            except:
                heights_ = [heights]
        else:
            heights_ = self._heights

        for height in heights_:
            if fitreport or verbose:
                print('PSD measurment at {0:1.3f} {1:s}'
                      ''.format(height, self._height_unit))

            self.psdfits[height].fit_psds(names=names,
                                          dynamic_bounds=dynamic_bounds,
                                          bounds=bounds,
                                          f_exclude=f_exclude,
                                          use_heights=use_heights,
                                          calc_results=True,
                                          fitreport=fitreport,
                                          verbose=verbose,
                                          **kwargs)
            if plot_fits:
                plot_axis.cla()
                ttl = ('PSDs at {0:1.3f} {1:s}'
                       ''.format(height, self._height_unit))
                if save_plots:
                    suffix = int(height * conv)
                    fname = '{0:s}-h_{1:+d}nm'.format(bn, suffix)
                else:
                    fname = None

                fig = self.psdfits[height].plot_fits(names=names,
                                                     plot_axis=plot_axis,
                                                     plot_data=True,
                                                     plot_masked_data=True,
                                                     save_plots=save_plots,
                                                     save_as=save_as,
                                                     directory=directory,
                                                     filename=fname,
                                                     title=ttl,
                                                     **plot_kws)
                fig.canvas.draw()
                time.sleep(plot_delay)

        # check if all psdfits were fitted
        if all(psdfit.pc_results != {} for psdfit in self.psdfits.values()):
            self._gen_height_data()
            self._construct_arrays()

        if plot_fits:
            plt.ioff()

    #####----------------------------------------------------------------------
    #---- height dependent data ----
    #####----------------------------------------------------------------------

    def _gen_height_data(self):
        """
        Generate the data vectors according to the heights.
        """
        self.sort_by_height()

        h_ = list(self.get_heights(unit=self._height_unit))
        pf = self.psdfits

        dissens_unit = self._dissens_unit
        k_unit = self._trap_stiffness_unit
        drag_unit = self._drag_unit

        ex_axis = 'excited_axis'

        for name in self._names:
            key = 'beta_' + name + '_pc'
            val = array([pf[h].pc_results[name].get_dissens(unit=dissens_unit)
                         for h in h_])
            self._dissens[key] = val

            key = 'dbeta_' + name + '_pc'
            val = array([pf[h].pc_results[name].get_dissens_err(unit=dissens_unit)
                         for h in h_])
            self._dissens[key] = val

            key = 'beta_' + name + '_ac'
            val = array([pf[h].ac_results[name].get_dissens(unit=dissens_unit)
                         for h in h_])
            self._dissens[key] = val

            key = 'dbeta_' + name + '_ac'
            val = array([pf[h].ac_results[name].get_dissens_err(unit=dissens_unit)
                         for h in h_])
            self._dissens[key] = val

            key = 'kappa_' + name + '_pc'
            val = array([pf[h].pc_results[name].get_trap_stiffness(unit=k_unit)
                         for h in h_])
            self._trap_stiffness[key] = val

            key = 'dkappa_' + name + '_pc'
            val = array([pf[h].pc_results[name].get_trap_stiffness_err(unit=k_unit)
                         for h in h_])
            self._trap_stiffness[key] = val

            key = 'kappa_' + name + '_ac'
            val = array([pf[h].ac_results[name].get_trap_stiffness(unit=k_unit)
                         for h in h_])
            self._trap_stiffness[key] = val
            key = 'dkappa_' + name + '_ac'
            val = array([pf[h].ac_results[name].get_trap_stiffness_err(unit=k_unit)
                         for h in h_])
            self._trap_stiffness[key] = val

            key = 'red_chi2_' + name
            self._red_chi2[key] = array([pf[h].fits[name].redchi2 for h in h_])

            key = 'N_free_' + name
            self._red_chi2[key] = array([pf[h].fits[name].nfree for h in h_])

            key = 'offset_' + name
            self._offset[key] = array([pf[h].psdm.psds[name].offset for h in h_])

        self.generate_outliers_masks(conf_level=self.conf_level)

        dgs = array([pf[h].pc_results[ex_axis].get_drag(unit=drag_unit)
                     for h in h_])
        self._drag['drag_stokes'] = dgs

        ddgs = array([pf[h].pc_results[ex_axis].get_drag_err(unit=drag_unit)
                     for h in h_])
        self._drag['ddrag_stokes'] = ddgs

        dg = array([pf[h].ac_results[ex_axis].get_drag(unit=drag_unit)
                    for h in h_])
        self._drag['drag'] = dg

        ddg = array([pf[h].ac_results[ex_axis].get_drag_err(unit=drag_unit)
                     for h in h_])
        self._drag['ddrag'] = ddg

        rdg = array(dg / dgs)
        self._drag['rel_drag'] = rdg
        self._drag['drel_drag'] = array(rdg * ddg / dg)

        self._temp = array([pf[h].psdm.exp_setting.get_temp(self._temp_unit) for h in h_])
        self._temp_err = array([pf[h].psdm.exp_setting.get_temp_err(self._temp_unit)
                                for h in h_])
        self._ex_axes = array([self._names.index(pf[h].psdm.ex_axis)
                               for h in h_])

    def get_dissens(self,
                    name,
                    ac_data=True,
                    error=False,
                    unit=None,
                    getQuantity=False,
                    nomask=False,
                    get_outliers=False,
                    inverse_mask=False
                    ):
        """
        Return the height dependent displacement sensitivity data or their
        errors according to ac or pc calculations.

        Arguments
        ---------
        name : str
            Defines the name of the axis to be returned.
        ac_data : bool
            Return data for according to ac (True) or pc (False) calculations.
        error : bool
            Return the error instead.
        unit : str
            Defines the unit of the data. If None *self._dissens_unit* is used.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name=name,
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)

        if error:
            key = 'dbeta_'
        else:
            key = 'beta_'

        if ac_data:
            key += name + '_ac'
        else:
            key += name + '_pc'

        if unit is None:
            unit = self._dissens_unit
            output = self._dissens[key][~mask]
        else:
            conv = ureg(self._dissens_unit).to(unit).magnitude
            output = self._dissens[key][~mask] * conv

        if getQuantity:
            return output * ureg(unit)
        else:
            return output

    def set_dissens_unit(self, unit, convert=True):
        """
        Set the unit of the displacement sensitivity vectors. If convert is
        True the values are converted into this unit, if False their magnitude
        is not changed, but the unit.
        """
        unit = str(ureg(unit).units)
        for k, v in self._dissens.items():
            if convert:
                conv = ureg(self._dissens_unit).to(unit).magnitude
                self._dissens[k] = v * conv
        self._dissens_unit = unit

    def get_trap_stiffness(self,
                           name,
                           ac_data=True,
                           error=False,
                           unit=None,
                           getQuantity=False,
                           nomask=False,
                           get_outliers=False,
                           inverse_mask=False
                           ):
        """
        Return the height dependent trap_stiffness data or their
        errors according to ac or pc calculations.

        Arguments
        ---------
        name : str
            Defines the name of the axis to be returned.
        ac_data : bool
            Return data for according to ac (True) or pc (False) calculations.
        error : bool
            Return the error instead.
        unit : str
            Defines the unit of the data. If None *self._trap_stiffness_unit*
            is used.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name=name,
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)

        if error:
            key = 'dkappa_'
        else:
            key = 'kappa_'
        if ac_data:
            key += name + '_ac'
        else:
            key += name + '_pc'
        if unit is None:
            unit = self._trap_stiffness_unit
            output = self._trap_stiffness[key][~mask]
        else:
            conv = ureg(self._trap_stiffness_unit).to(unit).magnitude
            output = self._trap_stiffness[key][~mask] * conv
        if getQuantity:
            return output * ureg(unit)
        else:
            return output

    def set_trap_stiffness_unit(self, unit, convert=True):
        """
        Set the unit of the trap stiffness vectors. If convert is True
        the values are converted into this unit, if False their magnitude is
        not changed, but the unit.
        """
        unit = str(ureg(unit).units)
        for k, v in self._trap_stiffness.items():
            if convert:
                conv = ureg(self._trap_stiffness_unit).to(unit).magnitude
                self._trap_stiffness[k] = v * conv
        self._trap_stiffness_unit = unit

    def get_drag(self,
                 error=False,
                 relative_drag=False,
                 stokes=False,
                 unit=None,
                 getQuantity=False,
                 nomask=False,
                 get_outliers=False,
                 inverse_mask=False):
        """
        Return the height dependent drag, stokes drag or relative drag data or
        their errors.

        Arguments
        ---------
        relative_drag : bool
            Return the relative drag if True.
        stokes : bool
            Return the stokes drag if True, the measured drag if False.
        error : bool
            Return the error instead.
        unit : str
            Defines the unit of the data. If None *self._drag_unit*
            is used.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name='ex_axis',
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)
        if error:
            key = 'd'
        else:
            key = ''

        if relative_drag:
            key += 'rel_drag'
            unit = 'dimensionless'
            if getQuantity:
                return self._drag[key][~mask] * ureg(unit)
            else:
                return self._drag[key][~mask]

        if stokes:
            key += 'drag_stokes'
        else:
            key += 'drag'

        if unit is None:
            unit = self._drag_unit
            output = self._drag[key][~mask]
        else:
            conv = ureg(self._drag_unit).to(unit).magnitude
            output = self._drag[key][~mask] * conv
        if getQuantity:
            return output * ureg(unit)
        else:
            return output

    def set_drag_unit(self, unit, convert=True):
        """
        Set the unit of the drag vectors. If convert is True
        the values are converted into this unit, if False their magnitude is
        not changed, but the unit.
        """
        unit = str(ureg(unit).units)
        for k, v in self._drag.items():
            if k.find('rel_') < 0:
                if convert:
                    conv = ureg(self._drag_unit).to(unit).magnitude
                    self._drag[k] = v * conv
        self._drag_unit = unit

    def get_redchi2(self,
                    name,
                    nomask=False,
                    get_outliers=False,
                    inverse_mask=False
                    ):
        """
        Return the height dependent displacement sensitivity data or their
        errors according to ac or pc calculations.

        Arguments
        ---------
        name : str
            Defines the name of the axis to be returned.
        ac_data : bool
            Return data for according to ac (True) or pc (False) calculations.
        error : bool
            Return the error instead.
        unit : str
            Defines the unit of the data. If None *self._dissens_unit* is used.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name=name,
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)

        return self._red_chi2['red_chi2_' + name][~mask]

    def get_offset(self,
                   name,
                   nomask=False,
                   get_outliers=False,
                   inverse_mask=False
                   ):
        """
        Return the height-dependent offset (dc-psd) signal in units of the psd.

        Arguments
        ---------
        name : str
            Defines the name of the axis to be returned.
        nomask : bool
            Whether to return a masked array or the full data vector.
        """
        if nomask:
            mask = zeros(shape(self._heights)) > 0
        else:
            mask = self.get_mask(name=name,
                                 get_outliers=get_outliers,
                                 inverse=inverse_mask)

        return self._offset['offset_' + name][~mask]

    #####----------------------------------------------------------------------
    #---- save and load height dependent data ----
    #####----------------------------------------------------------------------

    def _construct_arrays(self, get_all=False):
        """
        Gathers the values from the *PSDFit* objects into arrays.

        array names for names = ['x', 'y', 'z']
         - height
         - drag_stokes ...
         - ddrag_stokes ...
         - beta_x_pc, beta_y_pc, beta_z_pc
         - dbeta_x_pc ...
         - kappa_x_pc ...
         - dkappa_x_pc...
         - red_chi2_x ...
         - N_free_x ...
         - outlier_x ...
         - beta_x_ac ...
         - dbeta_x_ac ...
         - kappa_x_ac ...
         - dkappa_x_ac...
         - drag
         - ddrag
         - rel_drag
         - drel_drag
         - temp
        if get_all is True these arrays are added as well:
         - D_x ...
         - dD_x ...
         - f_c_x ...
         - df_c_x ...
         - A and dA - excitation amplitude and error
         - x- y- and z-signal (offsets)
         - radius
         - radius_err
        """
        a = OrderedDict()
        a[co.height] = self.get_heights(unit=self._height_unit, nomask=True)
        a.update(self._dissens)
        a.update(self._trap_stiffness)
        a.update(self._drag)
        a.update(self._red_chi2)
        a.update(self._offset)
        for name in self._names:
            a['outlier_' + name] = self._outlier_mask[name]
        a['temp'] = self._temp
        a['dtemp'] = self._temp_err
        a['ex_axis'] = self._ex_axes

        if get_all and len(self.psdfits) > 0:
            pf = self.psdfits
            h_ = list(self.get_heights(unit=self._height_unit))
            for name in self._names:
                a['D_' + name] = array([pf[h].fits[name].D for h in h_])
                a['dD_' + name] = array([pf[h].fits[name].D_err for h in h_])
                a['f_c_' + name] = array([pf[h].fits[name].f_c for h in h_])
                a['df_c_' + name] = array([pf[h].fits[name].f_c_err
                                          for h in h_])
            a['ex_amp'] = array([pf[h].psdm.ex_amplitude for h in h_])
            a['dex_amp'] = array([pf[h].psdm.ex_amplitude_err for h in h_])
            a['P_ex'] = array([pf[h].psdm.ex_power for h in h_])
            a['dP_ex'] = array([pf[h].psdm.ex_power_err for h in h_])
            a['P_base'] = array([pf[h].get_basepower()[0] for h in h_])
            a['dP_base'] = array([pf[h].get_basepower()[1] for h in h_])

        self._arrays = a

    def save_arrays_to_file(self,
                            basename,
                            directory=None,
                            suffix='_hc_data',
                            extension='.dat',
                            append_only=False,
                            report=True
                            ):
        """
        Save the height-dependent values given in the *_arrays* dictionary to a
        file *basename_hc_data.dat*.
        """
        if directory is None:
            try:
                directory = self.directory
            except:
                directory = './'

        path = join(directory, basename + suffix + extension)

        data = array([self._arrays[k]
                      for k in self._arrays.keys()]).transpose()

        if append_only:
            write_opt = 'a'
        else:
            write_opt = 'w'

        with open(path, write_opt) as fl:
            #fl = open(path, 'w')
            fl.write('\t'.join(self._arrays.keys()) + '\n')
            for i in range(len(data)):
                s = ['{0:1.5E}'.format(row) for row in data[i]]
                fl.write('\t'.join(s) + '\n')

        if report:
            print('Data written to\t\t {}'.format(path))

    def load_arrays(self, path=None):
        """
        Load the height-dependent-data from the given file to the *arrays*
        dictionary.
        """
        if path is None:
            try:
                path = join(self.directory, self.datafilename)
            except:
                raise Exception('File not found at {}'
                                ''.format(join(self.directory,
                                               self.datafilename)))

        arrays = read_std_data_file(path)
        self._arrays = arrays
        self._heights = list(arrays[co.height])
        for k, val in arrays.items():
            if k.startswith('beta_') or k.startswith('dbeta_'):
                self._dissens[k] = val
            elif k.startswith('kappa_') or k.startswith('dkappa_'):
                self._trap_stiffness[k] = val
            elif k.find('drag') >= 0:
                self._drag[k] = val
            elif k.startswith('red_chi2_') or k.startswith('N_free_'):
                self._red_chi2[k] = val
            elif k.startswith('offset_'):
                self._offset[k] = val
            elif k == 'temp':
                self._temp = val
            elif k == 'dtemp':
                self._temp_err = val
            elif k == 'ex_axis':
                self._ex_axes = val

        self.reset_mask()
        self.generate_outliers_masks()

    def save_parameter_file(self,
                            basename,
                            directory=None,
                            suffix='_hc_parameters',
                            extension='.txt',
                            append_only=False,
                            section='DEFAULT',
                            report=True
                            ):
        """
        Save the parameters to a parameter file.
        """
        if directory is None:
            try:
                directory = self.directory
            except:
                directory = './'

        ptf = join(directory, basename + suffix + extension)

        d = OrderedDict()
        d[co.names] = ','.join(self.names)
        d[co.laterality] = ','.join(str(1 * self._laterality[name])
                                    for name in self.names)
        d[co.height_unit] = self._height_unit
        d['specified_radius'] = self.get_radius(unit=self._radius_unit)
        d['specified_radius_error'] = self.get_radius_err(unit=self._radius_unit)
        d[co.radius_unit] = self._radius_unit
        d[co.height_offset] = self.height_offset
        d[co.wl] = self.get_wavelength(unit='nm')
        d[co.wl_unit] = 'nm'
        d[co.fs] = self.focal_shift
        d[co.dfs] = self.focal_shift_err
        d[co.ref_ind] = self.ref_ind
        d[co.conf_level] = self.conf_level
        d[co.beta_unit] = self._dissens_unit
        d[co.kappa_unit] = self._trap_stiffness_unit
        d[co.drag_unit] = self._drag_unit
        d[co.temp_unit] = self._temp_unit

        pars = configparser.ConfigParser()
        pars[section] = d

        if append_only:
            write_opt = 'a'
        else:
            write_opt = 'w'

        with open(ptf, write_opt) as pfile:
            pars.write(pfile)

        if report:
            print('Parameters written to\t {}'.format(ptf))

    def load_parameter_file(self, path=None, section='DEFAULT'):
        """
        Load parameters from the given parameter file.
        """
        if path is None:
            path = join(self.directory, self.paramfilename)

        pars = configparser.ConfigParser()

        if not pars.read(path):
            raise Exception('Could not read parameter file {0:s}.'
                            ''.format(path))
        d = pars[section]

        names = d[co.names].replace(' ', '').split(',')
        self._names = sorted(names)

        lats = d[co.laterality].replace(' ', '').split(',')
        self._laterality = {names[i]: (lat == '1') for i, lat in enumerate(lats)}
        self._height_unit = str(ureg(d[co.height_unit]).units)
        self.set_radius(float(d['specified_radius']), unit=d[co.radius_unit])
        self.set_radius_err(float(d['specified_radius_error']),
                            unit=d[co.radius_unit])
        self.height_offset = float(d[co.height_offset])
        self.set_wavelength(float(d[co.wl]), unit=d[co.wl_unit])
        self.focal_shift = float(d[co.fs])
        self.focal_shift_err = float(d[co.dfs])
        self.ref_ind = float(d[co.ref_ind])
        self.conf_level = float(d[co.conf_level])
        self._dissens_unit = d[co.beta_unit]
        self._trap_stiffness_unit = d[co.kappa_unit]
        self._drag_unit = d[co.drag_unit]
        self._temp_unit = d[co.temp_unit]

    def save_hc_data(self,
                     basename=None,
                     directory=None,
                     dfile_suffix='_hc_data',
                     pfile_suffix='_hc_parameters',
                     get_all=False,
                     report=True
                     ):
        """
        Save the height calibration data to text files.

        The Result from the psd fits and the parameters are saved to
        text files.
        """
        if get_all and len(self.psdfits) == 0:
            warnings.warn('The "get_all" option is only available if all PSD '
                          'measurements were loaded and fitted before.')
            get_all = False
        self._construct_arrays(get_all=get_all)

        if basename is None:
            if self.basename:
                basename = self.basename
            else:
                basename = 'height_calibration'

        self.save_arrays_to_file(basename,
                                 directory=directory,
                                 suffix=dfile_suffix)
        self.save_parameter_file(basename,
                                 directory=directory,
                                 suffix=pfile_suffix)

    def load_hc_data(self,
                     basename,
                     directory,
                     dfile_suffix='_hc_data',
                     pfile_suffix='_hc_parameters',
                     section='DEFAULT'):
        self.basename = basename
        self.directory = directory
        self.paramfilename = basename + pfile_suffix + '.txt'
        self.datafilename = basename + dfile_suffix + '.dat'

        self.load_parameter_file(path=join(directory, self.paramfilename))
        self.load_arrays(path=join(directory, self.datafilename))

    #####----------------------------------------------------------------------
    #---- fit height dependent data ----
    #####----------------------------------------------------------------------

    def determine_focal_shift(self,
                              name='z',
                              signal=None,
                              idx_slice=None,
                              wavelength=None,
                              ref_index=1.326,
                              report_fs=True,
                              plot_fit=True,
                              wavelength_unit='nm',
                              height_unit=None,
                              **kws):
        """
        Try to determine the focal shift from the oscillation of the signal.

        Uses pyotc.focal_shift.get_focal_shift(). 

        Arguments
        ---------
        name : str
            Name of the axis with the oscillations.
        signal : array or None
            The signal with the oscillations with same length as heights of the
            given axis name. if None, the displacement sensitivity of the given
            axis name is used.
        idx_slice : slice
            Slice that defines the used data.
        wavelength : float
            Wavelength of the trapping laser. use the keyword argument
            wavelength_unit to define the unit. defualt is 'nm'
        ref_index : float
            Refractive index of the medium.
        report_fs : bool
            Wheter to print the found focal shift.
        plot_fit : bool

        Notes 
        -----
        The used data does not take the unmasked data, which is defined through
        HC._mask. Instead the idx_slice can be used to define the used data.
        """
        heights = self.get_heights(name=name, nomask=True,
                                   unit=height_unit or self._height_unit)
        if not idx_slice:
            idx_slice = slice(len(heights))

        if signal is None:
            signal = self.get_dissens(name, nomask=True)

        wl = wavelength or self.get_wavelength(unit=wavelength_unit)
        
        mzr = get_focal_shift(heights[idx_slice],
                              signal[idx_slice],
                              wl,
                              ref_index,
                              report_focal_shift=report_fs,
                              plot_fit=plot_fit,
                              wavelength_unit=wavelength_unit,
                              height_unit=height_unit or self._height_unit,
                              **kws)
        return mzr

    def fit_rel_drag(self,
                     method='radius',
                     lateral=True,
                     parameters=None,
                     plot_fit=False,
                     verbose=False,
                     fit_focal_shift=False,
                     **init_kwargs):
        """
        Fit the relative drag with respect to height.

        Arguments
        ---------
        method : {'viscosity', 'radius'}
            Whether correcting for viscosity or radius.
            N.B. One can use the method 'viscosity' still to correct
            for deviations of the radius, if the deviations are small and one
            can be sure that the viscosity did not deviate, e.g. due to local
            heating or other ingredients in the solution.
        lateral : bool
            use True if excitation was parallel to surface.
        parameters : lmfit.Parameters
            Parameters object. See also gen_height_fit_pars().
        plot_fit : bool
        verbose : bool
            Print fit results.

        Keyword Arguments
        -----------------
        Keyword arguments defining the initial values of the fit, These are
        passed to gen_height_fit_pars(), e.g.
        h0 = 0.0, corr = 1.0, h0min = -1.0, etc.
        """
        heights = self.get_heights(name='ex_axis', unit='um')
        rel_drag = self.get_drag(relative_drag=True)
        errors = self.get_drag(relative_drag=True, error=True)
        radius = self.get_radius(unit='um')
        if 'focal_shift' not in init_kwargs:
            init_kwargs['focal_shift'] = self.focal_shift

        drag_minimizer = fit_rel_drag(heights,
                                      rel_drag,
                                      errors,
                                      radius,
                                      parameters=parameters,
                                      lateral=lateral,
                                      method=method,
                                      fit_focal_shift=fit_focal_shift,
                                      verbose=verbose,
                                      **init_kwargs)

        pars = drag_minimizer.params
        drag_fit_eval = _rel_drag_fit_fun(pars,
                                          heights,
                                          radius,
                                          fit_data=None,
                                          errors=None,
                                          lateral=lateral,
                                          method=method)
        drag_fit_resi = _rel_drag_fit_fun(pars,
                                          heights,
                                          radius,
                                          fit_data=rel_drag,
                                          errors=errors,
                                          lateral=lateral,
                                          method=method)
        drag_minimizer.eval = drag_fit_eval
        drag_minimizer.residuals = drag_fit_resi
        self.rel_drag_fit = drag_minimizer

        if plot_fit:
            fig = self.plot_drag()
            self.plot_rel_drag_fit(axis=fig.axes[0])
            return fig

    def _gen_pc_dicts(self):
        """
        Returns
        -------
        4-tuple of the following dictionaries:
        heights_dict : dict
        data_dict : dict
        errors_dict : dict
        lateral_dict : dict
            Dictionary with names of axes as key names and value = True, if
            this axis is parallel to the surface and False otherwise. Any 'z'
            within the name of the current axis will result in a False,
            meaning an axial flag.
        """
        heights_dict = OrderedDict()
        data_dict = OrderedDict()
        errors_dict = OrderedDict()
        lateral_dict = OrderedDict()

        heights_dict['rel_drag'] = self.get_heights(name='ex_axis')
        data_dict['rel_drag'] = self.get_drag(relative_drag=True)
        errors_dict['rel_drag'] = self.get_drag(relative_drag=True, error=True)

        for name in self._names:
            height_ = self.get_heights(name=name, unit='um')
            beta_ = self.get_dissens(name, ac_data=False,
                                     error=False, unit='nm/mV')
            kappa_ = self.get_trap_stiffness(name, ac_data=False,
                                             error=False, unit='pN/nm')
            dbeta_ = self.get_dissens(name, ac_data=False,
                                      error=True, unit='nm/mV')
            dkappa_ = self.get_trap_stiffness(name, ac_data=False,
                                              error=True, unit='pN/nm')

            heights_dict['beta_' + name + '_pc'] = height_
            heights_dict['kappa_' + name + '_pc'] = height_
            data_dict['beta_' + name + '_pc'] = beta_
            data_dict['kappa_' + name + '_pc'] = kappa_
            errors_dict['beta_' + name + '_pc'] = dbeta_
            errors_dict['kappa_' + name + '_pc'] = dkappa_
        lateral_dict = self._laterality
        lateral_dict['rel_drag'] = True

        return (heights_dict, data_dict, errors_dict, lateral_dict)

    def fit_height_data(self,
                        method='radius',
                        fit_drag_first=True,
                        fit_focal_shift=False,
                        fit_dissens_osci=False,
                        plot_fit=False,
                        parameters=None,
                        fit_report=True,
                        verbose=False,
                        **init_kwargs):
        """
        Fit the height-dependent passive calibration results.

        Arguments
        ---------
        method : {'viscosity', 'radius'}
            whether correcting for viscosity or radius.
            N.B. One can use the method 'viscosity' still to correct
            for deviations of the radius, if the deviations are small and one
            can be sure that the viscosity did not deviate, e.g. due to local
            heating or other ingredients in the solution.
        fit_drag_first : bool
            If True, fits the relative drag data first to find better initial
            values for the global fit.
        fit_dissens_osci : bool
            Whether to try to fit oscillations on the displacement sensitivity.
            If True, the fit will have an additional summand of the shape:
            A * exp(d * heights) *
            sin( 4 * pi / (n * wavelength) * focal_shift * heigths + phi).
        plot_fit : bool
            Plot the results.
        parameters : lmfit.Parameters()
            Parameters object holding the parameters of the height fit. These
            can be generated with gen_height_fit_pars().
        verbose : bool
            Print some info to stdio.

        Keyword Arguments
        -----------------
        Keyword arguments defining the initial values of the fit, These are
        passed to gen_height_fit_pars(), e.g.
        h0 = 0.0, corr = 1.0, h0min = -1.0, etc.
        """
        names = self._names
        radius = self.get_radius(unit='um')

        if 'focal_shift' not in init_kwargs:
            init_kwargs['focal_shift'] = self.focal_shift

        (heights_dict, data_dict,
         errors_dict, lateral_dict) = self._gen_pc_dicts()

        ## check if there's enough data
        # get shortest height vector
        ndata = min(len(h) for h in heights_dict.values())
        # calculate needed number = 4 * naxes + 2 (2x slope, 2x offset +
        # surface and correction factor)
        needed = len(self._names) * 4 + 2
        if ndata < needed:
            raise Exception('Not enough data points! {} are given, '
                            'at least {} are needed.\n'
                            'You can try to exclude a hole axis, if fit of a '
                            'particular have too many outliers.'
                            ''.format(ndata, needed))

        if fit_drag_first:
            self.fit_rel_drag(method=method,
                              fit_focal_shift=False,
                              **init_kwargs)

        if self.rel_drag_fit:
            if 'h0' not in init_kwargs:
                init_kwargs['h0'] = self.rel_drag_fit.params['h0'].value
            if 'corr' not in init_kwargs:
                init_kwargs['corr'] = self.rel_drag_fit.params['corr'].value

        if fit_dissens_osci:
            fit_fun_kw = {'ref_ind': self.ref_ind,
                          'wavelength': self.get_wavelength(unit='um')}
        else:
            fit_fun_kw = {}

        hmin, dmin = fit_height_data(names,
                                     heights_dict,
                                     data_dict,
                                     errors_dict,
                                     radius,
                                     parameters=parameters,
                                     fit_drag_first=fit_drag_first,
                                     drag_fit_method=method,
                                     lateral_dict=lateral_dict,
                                     fit_dissens_osci=fit_dissens_osci,
                                     fit_fun_kw=fit_fun_kw,
                                     fit_focal_shift=fit_focal_shift,
                                     iterations=1,
                                     fit_report=fit_report,
                                     verbose=verbose,
                                     **init_kwargs)

        rel_drag_fun = (lambda h_, lateral:
                        _rel_drag_fit_fun(hmin.params,
                                          h_,
                                          radius,
                                          fit_data=None,
                                          errors=None,
                                          lateral=lateral,
                                          method=method))
        hmin.rel_drag_fun = rel_drag_fun

        hmin.eval = _height_fit_fun(hmin.params,
                                    names,
                                    heights_dict,
                                    radius,
                                    data_dict=None,
                                    errors_dict=None,
                                    lateral_dict=lateral_dict,
                                    method=method,
                                    return_fit_dict=True,
                                    fit_dissens_osci=fit_dissens_osci,
                                    dissens_fun_kw=fit_fun_kw)

        hmin.residuals = _height_fit_fun(hmin.params,
                                         names,
                                         heights_dict,
                                         radius,
                                         data_dict=data_dict,
                                         errors_dict=errors_dict,
                                         lateral_dict=lateral_dict,
                                         method=method,
                                         return_fit_dict=True,
                                         fit_dissens_osci=fit_dissens_osci,
                                         dissens_fun_kw=fit_fun_kw)
        self.height_fit = hmin

        hfr = HeightFitResult(method=method, minimizer=hmin)
        self.height_fit_results = hfr

        if plot_fit:
            plt.ion()
            self.plot_results()
            plt.ioff()

        return hfr

    def write_results_to_file(self,
                              basename='',
                              directory=None,
                              suffix='_hc_results',
                              extension='.txt',
                              section='DEFAULT',
                              append_only=False,
                              radius_unit='um',
                              height_unit='um',
                              dissens_unit='nm/mV',
                              dissens_slope_unit='nm/(mV*um)',
                              trap_stiffness_unit='pN/nm',
                              trap_stiffness_slope_unit='pN/(nm*um)'
                              ):
        """
        Write the results from the height fit to a file.
        """
        data = OrderedDict()
        pars = self.height_fit_results.params

        data[co.fit_method] = self.height_fit_results.method
        data[co.radius] = '{0:1.5e}'.format(self.get_radius(unit=radius_unit))
        data[co.radius_err] = '{0:1.5e}'.format(self.get_radius_err(unit=radius_unit))
        
        data[co.fs] = '{0:1.5e}'.format(pars['focal_shift'].value)
        data[co.dfs] = '{0:1.5e}'.format(pars['focal_shift'].stderr)

        if height_unit != 'um':
            hconv = ureg('um').to(height_unit).magnitude
        else:
            hconv = 1.0
        
        data[co.h0] = '{0:1.5e}'.format(pars['h0'].value * hconv)
        data[co.h0_err] = '{0:1.5e}'.format(pars['h0'].stderr * hconv)
        
        data[co.height_offset] = '{0:1.5e}'.format(self.height_offset * hconv)
        
        d_app_surf = (pars['h0'].value - self.height_offset) * hconv
        data[co.d_app_surf] = '{0:1.5e}'.format(d_app_surf)
        data[co.d_app_surf_err] = '{0:1.5e}'.format(pars['h0'].stderr * hconv)

        data[co.corr] = '{0:1.5e}'.format(pars['corr'].value)
        data[co.corr_err] = '{0:1.5e}'.format(pars['corr'].stderr)

        if dissens_unit != 'nm/mV':
            beta_conv = ureg('nm/mV').to(dissens_unit).magnitude
        else:
            beta_conv = 1.0

        if dissens_slope_unit != 'nm/(mV*um)':
            mbeta_conv = ureg('nm/(mV*um)').to(dissens_slope_unit).magnitude
        else:
            mbeta_conv = 1.0

        if trap_stiffness_unit != 'pN/nm':
            kappa_conv = ureg('pN/nm').to(trap_stiffness_unit).magnitude
        else:
            kappa_conv = 1.0

        if trap_stiffness_slope_unit != 'pN/(nm*um)':
            mkappa_conv = ureg('pN/(nm*um)').to(trap_stiffness_slope_unit).magnitude
        else:
            mkappa_conv = 1.0

        for name in self.names:
            data['beta_' + name] = '{0:1.5e}'.format(pars['beta_' + name].value * beta_conv)
            data['beta_' + name + '_err'] = '{0:1.5e}'.format(pars['beta_' + name].stderr * beta_conv)
            data['mbeta_' + name] = '{0:1.5e}'.format(pars['mbeta_' + name].value * mbeta_conv)
            data['mbeta_' + name + '_err'] = '{0:1.5e}'.format(pars['mbeta_' + name].stderr * mbeta_conv)

            data['kappa_' + name] = '{0:1.5e}'.format(pars['kappa_' + name].value * kappa_conv)
            data['kappa_' + name + '_err'] = '{0:1.5e}'.format(pars['kappa_' + name].stderr * kappa_conv)
            data['mkappa_' + name] = '{0:1.5e}'.format(pars['mkappa_' + name].value * mkappa_conv)
            data['mkappa_' + name + '_err'] = '{0:1.5e}'.format(pars['mkappa_' + name].stderr * mkappa_conv)

        units = OrderedDict()
        units[co.radius] = radius_unit
        units['apparent_distance'] = height_unit
        units['beta'] = dissens_unit
        units['mbeta'] = dissens_slope_unit
        units['kappa'] = trap_stiffness_unit
        units['mkappa'] = trap_stiffness_slope_unit

        parameters = configparser.ConfigParser()
        parameters['RESULTS'] = data
        parameters['UNITS'] = units

        if directory is None:
            try:
                directory = self.directory
            except:
                directory = './'
        if basename is '':
            try:
                basename = self.basename
            except:
                basename = ''

        if append_only:
            write_opt = 'a'
        else:
            write_opt = 'w'

        ptf = join(directory, basename + suffix + extension)

        with open(ptf, write_opt) as pfile:
            parameters.write(pfile)

        return ptf

    #####----------------------------------------------------------------------
    #---- plotting ----
    #####----------------------------------------------------------------------

    def plot_a_vs_b(self, name1, name2=co.height, fig=None, axis=None,
                    xlabel=None, ylabel=None, **kwargs):
        """
        Plot values of *name1* vs the values of *name2*.
        """
        if len(self._arrays) == 0:
            self._construct_arrays()
        a = self._arrays[name1]
        b = self._arrays[name2]
        if xlabel is None:
            xlabel = name2
        if ylabel is None:
            ylabel = name1
        ax = add_plot_to_figure(fig, b, a,
                                axis=axis,
                                xlabel=xlabel,
                                ylabel=ylabel,
                                **kwargs
                                )
        return ax.figure

    def plot_drag(self,
                  relative=True,
                  plot_errors=True,
                  axis=None,
                  nomask=False,
                  plot_outliers=False,
                  inverse_mask=False,
                  height_unit=None,
                  drag_unit=None,
                  **plot_kwargs):
        """
        Plot the measured drag vs. height.

        Arguments
        ---------
        relative : bool
            Plot the relative drag if True.
        plot_errors : bool
            Plot the error bars.
        axis : Axis
            matplotlib.Axis to add the plot to.
        nomask : bool
            plot all data if True.
        plot_outliers : bool
            Plot data points outside the confidence interval.
        inverse_mask : bool
            Whether to inverse the current mask.
        height_unit : str or None
            Unit of the height vector.
        drag_unit : str or None
            Unit of the drag (if relative is False). If None, self._drag_unit
            is used.
        """
        ax = axis
        if axis is None:
            fig = None
        else:
            fig = axis.figure

        if height_unit is None:
            height_unit = self._height_unit
        if drag_unit is None:
            drag_unit = self._drag_unit

        if 'color' not in plot_kwargs:
            if plot_outliers:
                plot_kwargs['color'] = col_dict['odrag']
            else:
                plot_kwargs['color'] = col_dict['drag']

        heights = self.get_heights(name='ex_axis',
                                   unit=height_unit,
                                   nomask=nomask,
                                   get_outliers=plot_outliers,
                                   inverse_mask=inverse_mask)

        if plot_errors:
            errors = self.get_drag(relative_drag=True,
                                   error=True,
                                   unit=drag_unit,
                                   nomask=nomask,
                                   get_outliers=plot_outliers,
                                   inverse_mask=inverse_mask)
        else:
            errors = []

        drag = self.get_drag(relative_drag=relative,
                             nomask=nomask,
                             get_outliers=plot_outliers,
                             inverse_mask=inverse_mask)

        if height_unit is None:
            height_unit = self._height_unit
        if drag_unit is None:
            drag_unit = self._drag_unit

        label = 'Relative drag'
        if plot_outliers:
            label += ' outliers'
        if inverse_mask:
            label += ' inverse'

        if relative:
            ylabel = '$\mathsf{{ \gamma / \gamma_0 }}$'
        else:
            ylabel = '$\mathsf{{ \gamma \; ({}) }}$'.format(drag_unit)

        ax = add_plot_to_figure(fig,
                                heights,
                                drag,
                                yerr=errors,
                                axis=ax,
                                fmt='o',
                                markersize=5,
                                label=label,
                                alpha=0.7,
                                xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                        ''.format(u2str(height_unit))),
                                ylabel=ylabel,
                                **plot_kwargs)
        return ax.figure

    def plot_rel_drag_fit(self,
                          height_fit=False,
                          plot_residuals=False,
                          axis=None,
                          height_unit=None,
                          **plot_kwargs):
        """
        Plot the evaluated fit or the residuals.

        Arguments
        ---------
        height_fit : bool
            If True, evaluation from (global) height fit is used.
        plot_residuals : bool
            Plot residuals instead of fit.
        plot_errors : bool
            Whether to plot the errors instead of the residuals.
        axis : matplotlib.Axis
            Axis to add the plot to.
        height_unit : str or None
            unit of the height vector.
        **plot_kwargs : keyword arguments
            passed over to plot function.
        """
        ax = axis
        if ax is None:
            fig = None
        else:
            fig = ax.figure

        heights = self.get_heights(name='ex_axis', unit=height_unit)

        if height_unit is None:
            height_unit = self._height_unit

        if 'color' not in plot_kwargs:
            plot_kwargs['color'] = col_dict['drag']

        if plot_residuals:
            if height_fit:
                ydata = self.height_fit.residuals['rel_drag']
                label = 'Residuals of global relative-drag-fit'
            else:
                ydata = self.rel_drag_fit.residuals
                label = 'Residuals of relative-drag-fit'

            '''
            if plot_errors:
                yerr = self.get_drag(relative_drag=True, error=True)
                ydata *= yerr
                ylabel = 'Errors'
            else:
                yerr = []
                ylabel='Residuals'
            '''

            ax = add_plot_to_figure(fig,
                                    heights,
                                    zeros(shape(heights)),
                                    fmt='-k',
                                    alpha=0.5,
                                    axis=ax)

            add_plot_to_figure(fig,
                               heights,
                               ydata,
                               axis=ax,
                               fmt='o',
                               markersize=5,
                               label=label,
                               alpha=0.7,
                               xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                       ''.format(u2str(height_unit))),
                               ylabel='Residuals',
                               **plot_kwargs)
            set_sym_y_labels(ax, mean=0)

            return ax.figure

        else:
            if height_fit:
                fit = self.height_fit.eval['rel_drag']
                label = 'Global relative-drag-fit'
            else:
                fit = self.rel_drag_fit.eval
                label = 'Relative-drag-fit'

            ax = add_plot_to_figure(fig,
                                    heights,
                                    fit,
                                    alpha=0.7,
                                    axis=ax,
                                    xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                            ''.format(u2str(height_unit))),
                                    ylabel='$\mathsf{{ \gamma / \gamma_0 }}$',
                                    label=label,
                                    **plot_kwargs)
            return ax.figure

    def plot_dissens(self,
                     plot_ac_data=True,
                     plot_errors=True,
                     plot_corrected=False,
                     names=None,
                     axis=None,
                     nomask=False,
                     plot_outliers=False,
                     inverse_mask=False,
                     height_unit=None,
                     dissens_unit=None,
                     **plot_kwargs):
        """
        Plot the measured displacement sensitivity vs. height.

        Arguments
        ---------
        plot_errors : bool
            Plot the error bars.
        name : str
            name(s) of the axis to be plotted or 'all' to plot all.
        axis : Axis
            matplotlib.Axis to add the plot to.
        nomask : bool
            plot all data if True.
        plot_outliers : bool
            Plot data points outside the confidence interval.
        inverse_mask : bool
            Whether to inverse the current mask.
        height_unit : str or None
            Unit of the height vector.
        dissens_unit : str
            Unit of the displacement sensitivity.
        """
        ax = axis
        if axis is None:
            fig = None
        else:
            fig = axis.figure

        if names and not isinstance(names, list):
            names = [names]

        if height_unit is None:
            height_unit = self._height_unit
        if dissens_unit is None:
            dissens_unit = self._dissens_unit

        if plot_ac_data and plot_corrected:
            raise Exception('There is no correction for the ac-data.')

        for idx, name in enumerate(names if names else self.names):
            heights = self.get_heights(name=name,
                                       unit=height_unit,
                                       nomask=nomask,
                                       get_outliers=plot_outliers,
                                       inverse_mask=inverse_mask)

            if plot_errors:
                errors = self.get_dissens(name=name,
                                          ac_data=plot_ac_data,
                                          unit=dissens_unit,
                                          error=True,
                                          nomask=nomask,
                                          get_outliers=plot_outliers,
                                          inverse_mask=inverse_mask)
            else:
                errors = []

            dissens = self.get_dissens(name=name,
                                       ac_data=plot_ac_data,
                                       unit=dissens_unit,
                                       error=False,
                                       nomask=nomask,
                                       get_outliers=plot_outliers,
                                       inverse_mask=inverse_mask)

            if plot_corrected:
                try:
                    lat = self._laterality[name]
                    rel_drag = self.height_fit.rel_drag_fun(heights, lat)
                except:
                    raise Exception('Before plotting the corrected data try'
                                    'to fit the height dependent data first.')
                data = dissens / rel_drag**0.5
            else:
                data = dissens

            if plot_ac_data:
                label = r'$\beta_{{\mathrm{{{ax},\, ac}}}}$'.format(ax=name)
            else:
                label = r'$\beta_{{\mathrm{{{ax},\, pc}}}}$'.format(ax=name)

            if plot_outliers:
                label += ' outliers'
            if inverse_mask:
                label += ' inverse'

            if name in col_dict.keys():
                if plot_outliers:
                    plot_kwargs['color'] = col_dict['o' + name]
                else:
                    plot_kwargs['color'] = col_dict[name]
            elif idx <= 3:
                if plot_outliers:
                    plot_kwargs['color'] = col_dict['o' + name]
                else:
                    plot_kwargs['color'] = col_dict[idx + 10]

            ax = add_plot_to_figure(fig,
                                    heights,
                                    data,
                                    yerr=errors,
                                    axis=ax,
                                    fmt='o',
                                    markersize=5,
                                    label=label,
                                    alpha=0.7,
                                    xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                            ''.format(u2str(height_unit))),
                                    ylabel=(r'$ \beta \; \mathrm{{({})}}$'
                                            ''.format(u2str(dissens_unit))),
                                    **plot_kwargs)
        return ax.figure

    def plot_dissens_fit(self,
                         names=None,
                         plot_residuals=False,
                         plot_corrected=False,
                         axis=None,
                         height_unit=None,
                         dissens_unit=None,
                         **plot_kwargs):
        """
        Plot the evaluated fit or the residuals.

        Arguments
        ---------
        height_fit : bool
            If True, evaluation from (global) height fit is used.
        plot_residuals : bool
            Plot residuals instead of fit.
        plot_errors : bool
            Whether to plot the errorbars.
        name : str
            name of the axis to be plotted or 'all' to plot all.
        axis : matplotlib.Axis
            Axis to add the plot to.
        height_unit : str or None
            unit of the height vector.
        **plot_kwargs : keyword arguments
            passed over to plot function.
        """
        ax = axis
        if ax is None:
            fig = None
        else:
            fig = ax.figure

        if height_unit is None:
            height_unit = self._height_unit
        if dissens_unit is None:
            dissens_unit = self._dissens_unit

        if names and not isinstance(names, list):
            names = [names]

        for idx, name in enumerate(names if names else self.names):
            heights = self.get_heights(name=name, unit=height_unit)

            if name in col_dict:
                plot_kwargs['color'] = col_dict[name]
            elif idx <= 3:
                plot_kwargs['color'] = col_dict[idx]

            if plot_residuals:
                ydata = self.height_fit.residuals['beta_' + name + '_pc']
                '''
                if plot_errors:
                    errors = self.get_dissens(name=name,
                                              error=True,
                                              ac_data=False)
                    ydata *= errors
                    ylabel='Errors'
                else:
                    errors = []
                    ylabel='Residuals'
                '''
                label = (r'$\beta_{{\mathrm{{{ax}}}}} \mathrm{{residuals}}$'
                         ''.format(ax=name))

                ax = add_plot_to_figure(fig,
                                        heights,
                                        zeros(shape(heights)),
                                        fmt='-k',
                                        alpha=0.5,
                                        axis=ax)

                add_plot_to_figure(fig,
                                   heights,
                                   ydata,
                                   axis=ax,
                                   fmt='o',
                                   markersize=5,
                                   label=label,
                                   alpha=0.7,
                                   xlabel=('Stage height ({})'
                                           ''.format(u2str(height_unit))),
                                   ylabel='Residuals',
                                   **plot_kwargs)
                set_sym_y_labels(ax, mean=0)
            else:
                fit = self.height_fit.eval['beta_' + name + '_pc']
                if plot_corrected:
                    try:
                        lat = self._laterality[name]
                        rel_drag = self.height_fit.rel_drag_fun(heights, lat)
                    except:
                        raise Exception('Before plotting the corrected data '
                                        'try to fit the height dependent data '
                                        'first.')
                    data = fit / rel_drag**0.5
                else:
                    data = fit

                if dissens_unit != 'nm/mV':
                    conv = ureg('nm/mV').to(dissens_unit).magnitude
                    data *= conv

                label = (r'$\beta_{{\mathrm{{{ax}, corrected}}}}$'
                         ''.format(ax=name))

                ax = add_plot_to_figure(fig,
                                        heights,
                                        data,
                                        axis=ax,
                                        fmt='-',
                                        label=label,
                                        alpha=0.7,
                                        xlabel=('Stage height ({})'
                                                ''.format(u2str(height_unit))),
                                        ylabel=(r'$ \beta \; \mathrm{{({})}}$'
                                                ''.format(u2str(dissens_unit))),
                                        **plot_kwargs)
        return ax.figure

    def plot_trap_stiffness(self,
                            plot_ac_data=True,
                            plot_errors=True,
                            plot_corrected=False,
                            names=None,
                            axis=None,
                            nomask=False,
                            plot_outliers=False,
                            inverse_mask=False,
                            height_unit=None,
                            trap_stiffness_unit=None,
                            **plot_kwargs):
        """
        Plot the measured trap stiffness vs. height.

        Arguments
        ---------
        plot_ac_data : bool
            Whether to plot pc or ac data.
        plot_errors : bool
            Plot the error bars.
        plot_corrected : bool
            Plot the pc measurements corrected by the drag.
        name : str
            name of the axis to be plotted or 'all' to plot all.
        axis : Axis
            matplotlib.Axis to add the plot to.
        nomask : bool
            plot all data if True.
        plot_outliers : bool
            Plot data points outside the confidence interval.
        inverse_mask : bool
            Whether to inverse the current mask.
        height_unit : str or None
            Unit of the height vector.
        trap_stiffness_unit : str
            Unit of the trap stiffness.
        """
        ax = axis
        if axis is None:
            fig = None
        else:
            fig = axis.figure

        if height_unit is None:
            height_unit = self._height_unit
        if trap_stiffness_unit is None:
            trap_stiffness_unit = self._trap_stiffness_unit

        if names and not isinstance(names, list):
            names = [names]

        if plot_ac_data and plot_corrected:
            raise Exception('There is no correction for the ac-data.')

        for idx, name in enumerate(names if names else self.names):
            heights = self.get_heights(name=name,
                                       unit=height_unit,
                                       nomask=nomask,
                                       get_outliers=plot_outliers,
                                       inverse_mask=inverse_mask)

            if plot_errors:
                errors = self.get_trap_stiffness(name=name,
                                                 ac_data=plot_ac_data,
                                                 unit=trap_stiffness_unit,
                                                 error=True,
                                                 nomask=nomask,
                                                 get_outliers=plot_outliers,
                                                 inverse_mask=inverse_mask)
            else:
                errors = []

            kappa = self.get_trap_stiffness(name=name,
                                            ac_data=plot_ac_data,
                                            unit=trap_stiffness_unit,
                                            error=False,
                                            nomask=nomask,
                                            get_outliers=plot_outliers,
                                            inverse_mask=inverse_mask)
            if plot_corrected:
                try:
                    lat = self._laterality[name]
                    rel_drag = self.height_fit.rel_drag_fun(heights, lat)
                except:
                    raise Exception('Before plotting the corrected data try'
                                    'to fit the height dependent data first.')
                data = kappa * rel_drag
            else:
                data = kappa

            if plot_ac_data:
                label = r'$\kappa_{{\mathrm{{{ax},\,ac}}}}$'.format(ax=name)
            else:
                label = r'$\kappa_{{\mathrm{{{ax},\,pc}}}}$'.format(ax=name)

            if plot_outliers:
                label += ' outliers'
            if inverse_mask:
                label += ' inverse'

            if name in col_dict.keys():
                if plot_outliers:
                    plot_kwargs['color'] = col_dict['o' + name]
                else:
                    plot_kwargs['color'] = col_dict[name]
            elif idx <= 3:
                if plot_outliers:
                    plot_kwargs['color'] = col_dict['o' + name]
                else:
                    plot_kwargs['color'] = col_dict[idx + 10]

            ax = add_plot_to_figure(fig,
                                    heights,
                                    data,
                                    yerr=errors,
                                    axis=ax,
                                    fmt='o',
                                    markersize=5,
                                    label=label,
                                    alpha=0.7,
                                    xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                            ''.format(u2str(height_unit))),
                                    ylabel=(r'$ \kappa \; \mathrm{{({})}}$'
                                            ''.format(u2str(trap_stiffness_unit))),
                                    **plot_kwargs)
        return ax.figure

    def plot_trap_stiffness_fit(self,
                                plot_residuals=False,
                                plot_corrected=False,
                                names=None,
                                axis=None,
                                height_unit=None,
                                trap_stiffness_unit=None,
                                **plot_kwargs):
        """
        Plot the evaluated fit or the residuals.

        Arguments
        ---------
        height_fit : bool
            If True, evaluation from (global) height fit is used.
        plot_residuals : bool
            Plot residuals instead of fit.
        plot_errors : bool
            Whether to plot the errorbars.
        name : str
            name of the axis to be plotted or 'all' to plot all.
        axis : matplotlib.Axis
            Axis to add the plot to.
        height_unit : str or None
            unit of the height vector.
        **plot_kwargs : keyword arguments
            passed over to plot function.
        """
        ax = axis
        if ax is None:
            fig = None
        else:
            fig = ax.figure

        if height_unit is None:
            height_unit = self._height_unit
        if trap_stiffness_unit is None:
            trap_stiffness_unit = self._trap_stiffness_unit

        if names and not isinstance(names, list):
            names = [names]

        for idx, name in enumerate(names if names else self.names):
            heights = self.get_heights(name=name, unit=height_unit)

            if name in col_dict:
                plot_kwargs['color'] = col_dict[name]
            elif idx <= 3:
                plot_kwargs['color'] = col_dict[idx]

            if plot_residuals:
                ydata = self.height_fit.residuals['kappa_' + name + '_pc']
                '''
                if plot_errors:
                    errors = self.get_trap_stiffness(name=name,
                                                     error=True,
                                                     ac_data=False)
                    ydata *= errors
                    ylabel = 'Errors'
                else:
                    errors = []
                    ylabel = 'Residuals'
                '''
                label = (r'$\kappa_{{\mathrm{{{ax}}}}} \mathrm{{residuals}}$'
                         ''.format(ax=name))

                ax = add_plot_to_figure(fig,
                                        heights,
                                        zeros(shape(heights)),
                                        fmt='-k',
                                        alpha=0.5,
                                        axis=ax)

                add_plot_to_figure(fig,
                                   heights,
                                   ydata,
                                   axis=ax,
                                   fmt='o',
                                   markersize=5,
                                   label=label,
                                   alpha=0.7,
                                   xlabel=('Stage height ({})'
                                           ''.format(u2str(height_unit))),
                                   ylabel='Residuals',
                                   **plot_kwargs)
                set_sym_y_labels(ax, mean=0)
            else:
                fit = self.height_fit.eval['kappa_' + name + '_pc']
                if plot_corrected:
                    try:
                        lat = self._laterality[name]
                        rel_drag = self.height_fit.rel_drag_fun(heights, lat)
                    except:
                        raise Exception('Before plotting the corrected data '
                                        'try to fit the height dependent data '
                                        'first.')
                    data = fit * rel_drag
                else:
                    data = fit

                if trap_stiffness_unit != 'pN/nm':
                    conv = ureg('pN/nm').to(trap_stiffness_unit).magnitude
                    data *= conv

                label = (r'$\kappa_{{\mathrm{{{ax}, corrected}}}}$'
                         ''.format(ax=name))

                add_plot_to_figure(fig,
                                   heights,
                                   data,
                                   axis=ax,
                                   fmt='-',
                                   label=label,
                                   alpha=0.7,
                                   xlabel=('Stage height ({})'
                                           ''.format(u2str(height_unit))),
                                   ylabel=(r'$ \kappa \; \mathrm{{({})}}$'
                                           ''.format(u2str(trap_stiffness_unit))),
                                   **plot_kwargs)
        return ax.figure

    def plot_redchi2(self,
                     names=None,
                     axis=None,
                     nomask=False,
                     plot_outliers=False,
                     inverse_mask=False,
                     height_unit=None,
                     **plot_kwargs):
        """
        Plot the measured displacement sensitivity vs. height.

        Arguments
        ---------
        plot_errors : bool
            Plot the error bars.
        name : str
            name(s) of the axis to be plotted or 'all' to plot all.
        axis : Axis
            matplotlib.Axis to add the plot to.
        nomask : bool
            plot all data if True.
        plot_outliers : bool
            Plot data points outside the confidence interval.
        inverse_mask : bool
            Whether to inverse the current mask.
        height_unit : str or None
            Unit of the height vector.
        dissens_unit : str
            Unit of the displacement sensitivity.
        """
        ax = axis
        if ax is None:
            fig = None
        else:
            fig = ax.figure

        if height_unit is None:
            height_unit = self._height_unit

        if names and not isinstance(names, list):
            names = [names]

        for idx, name in enumerate(names if names else self.names):
            heights = self.get_heights(name=name,
                                       unit=height_unit,
                                       nomask=nomask,
                                       get_outliers=plot_outliers,
                                       inverse_mask=inverse_mask)

            redchi2 = self.get_redchi2(name=name,
                                       nomask=nomask,
                                       get_outliers=plot_outliers,
                                       inverse_mask=inverse_mask)

            label = r'$\chi^2_{{\mathsf{{red,\; {ax} }} }}$'.format(ax=name)

            if plot_outliers:
                label += ' outliers'
            if inverse_mask:
                label += ' inverse'

            if 'c' not in plot_kwargs or 'color' not in plot_kwargs:
                if name in col_dict.keys():
                    plot_kwargs['color'] = col_dict[name]
                elif idx <= 3:
                    plot_kwargs['color'] = col_dict[idx]

            ax = add_plot_to_figure(fig,
                                    heights,
                                    redchi2,
                                    axis=ax,
                                    fmt='o',
                                    markersize=5,
                                    label=label,
                                    alpha=0.7,
                                    xlabel=('Stage height ($\mathsf{{ {} }}$)'
                                            ''.format(u2str(height_unit))),
                                    ylabel=(r'$ \chi^2_\mathsf{red}$'),
                                    **plot_kwargs)
        return ax.figure

    def plot_pc_results(self,
                        names=None,
                        plot_corrected=False,
                        figure=None,
                        figsize=(9, 9),
                        height_unit=None,
                        dissens_unit=None,
                        trap_stiffness_unit=None,
                        return_axes=False,
                        save=False,
                        path=None,
                        save_as='png',
                        dpi=150
                        ):
        """
        Plot the results from the fits to the power spectral densities.

        This function produces a graph with four diagrams showing the following
        values at a given height.
         - displacements sensitivities for x, y and z axis and the calibrated
           displacement sensitivity
         - trap stiffnesses for x, y and z and the trap stiffnes of the
           calibrated axis
         - reduced chi² values of the PSD fits
        """
        if names and not isinstance(names, list):
            names = [names]

        grid = (4, 1)

        if figure is None:
            fig = plt.figure(figsize=figsize)

        ax1 = get_gridded_axis(figure=fig,
                               grid=grid,
                               axslice=(slice(0, 1, 1), slice(0, 1, 1)))

        ax2 = get_gridded_axis(figure=fig,
                               grid=grid,
                               axslice=(slice(1, 2, 1),
                                        slice(0, 1, 1)),
                               sharex=ax1)

        ax3 = get_gridded_axis(figure=fig,
                               grid=grid,
                               axslice=(slice(2, 3, 1), slice(0, 1, 1)),
                               sharex=ax1)

        ax4 = get_gridded_axis(figure=fig,
                               grid=grid,
                               axslice=(slice(3, 4, 1), slice(0, 1, 1)),
                               sharex=ax1)

        ax2b = ax2.twinx()
        ax3b = ax3.twinx()
        ax4b = ax4.twinx()

        axes = {'drag': ax1,
                'dissens_left': ax2,
                'dissens_right': ax2b,
                'trap_stiffness_left': ax3,
                'trap_stiffness_right': ax3b,
                'red_chi2_left': ax4,
                'red_chi2_right': ax4b
                }

        self.plot_drag(axis=ax1, height_unit=height_unit)
        self.plot_drag(axis=ax1, plot_outliers=True,
                       height_unit=height_unit, c='gray')

        zplotted = False

        for idx, name in enumerate(names if names else self.names):
            if name.find('z') < 0:
                ax2_ = ax2
                ax3_ = ax3
                ax4_ = ax4
            else:
                zplotted = True
                ax2_ = ax2b
                ax3_ = ax3b
                ax4_ = ax4b

            self.plot_dissens(names=name,
                              axis=ax2_,
                              plot_ac_data=False,
                              plot_corrected=plot_corrected,
                              c=col_dict[idx],
                              height_unit=height_unit,
                              dissens_unit=dissens_unit)
            self.plot_dissens(names=name,
                              axis=ax2_,
                              plot_ac_data=False,
                              plot_corrected=plot_corrected,
                              plot_outliers=True,
                              c=col_dict[idx + 10],
                              height_unit=height_unit,
                              dissens_unit=dissens_unit)

            self.plot_trap_stiffness(names=name,
                                     axis=ax3_,
                                     plot_ac_data=False,
                                     plot_corrected=plot_corrected,
                                     c=col_dict[idx],
                                     height_unit=height_unit,
                                     trap_stiffness_unit=trap_stiffness_unit)
            self.plot_trap_stiffness(names=name,
                                     axis=ax3_,
                                     plot_ac_data=False,
                                     plot_corrected=plot_corrected,
                                     plot_outliers=True,
                                     c=col_dict[idx + 10],
                                     height_unit=height_unit,
                                     trap_stiffness_unit=trap_stiffness_unit)

            self.plot_redchi2(names=name, axis=ax4_, c=col_dict[idx])
            self.plot_redchi2(names=name, axis=ax4,
                              plot_outliers=True, c=col_dict[idx + 10])

        plt.setp([ax1, ax2, ax3], xlabel='')
        plt.setp([ax1.get_xticklabels(),
                  ax2.get_xticklabels(),
                  ax3.get_xticklabels()],
                 visible=False)

        if not zplotted:
            fig.delaxes(ax2b)
            fig.delaxes(ax3b)
            fig.delaxes(ax4b)

        fig.tight_layout()

        if save:
            if path is None:
                flnm = self.basename + '_results_plot.' + save_as
                path = join(self.directory, flnm)
            fig.patch.set_alpha(1.0)
            fig.patch.set_color('white')
            fig.savefig(path, dpi=dpi, format=save_as)

        if return_axes:
            return axes
        else:
            return fig

    def plot_rel_drag_fit_result(self,
                                 figure=None,
                                 figsize=(9, 6),
                                 save=False,
                                 path=None,
                                 save_as='png',
                                 dpi=150
                                 ):
        """
        Plot the relative drag, the fit and the residuals.
        """
        axes = get_residual_plot_axes(nrows=7,
                                      ncols=1,
                                      row_lim=5,
                                      figure=figure,
                                      figsize=figsize)
        ax1 = axes[0][0]
        ax2 = axes[0][1]
        fig = ax1.figure

        self.plot_drag(axis=ax1)
        self.plot_rel_drag_fit(axis=ax1)

        self.plot_rel_drag_fit(plot_residuals=True, axis=ax2)

        ax1.set_xlabel('')
        ax1.xaxis.set_ticklabels(ax1.get_xticklabels(), visible=False)

        fig.tight_layout()

        if save:
            if path is None:
                flnm = self.basename + '_results_plot.' + save_as
                path = join(self.directory, flnm)
            fig.patch.set_alpha(1.0)
            fig.patch.set_color('white')
            fig.savefig(path, dpi=dpi, format=save_as)

        return fig

    def plot_results(self,
                     names=None,
                     figure=None,
                     figsize=(14, 7),
                     plot_corrected=False,
                     height_unit=None,
                     dissens_unit=None,
                     trap_stiffness_unit=None,
                     save=False,
                     path=None,
                     save_as='png',
                     dpi=150):
        """
        DOC
        """
        if figure is None:
            figure = plt.figure(figsize=figsize)

        if names:
            if not isinstance(names, list):
                names = [names]
        else:
            names = self._names

        ax_drag, ax_dissens, ax_trap = get_residual_plot_axes(nrows=7,
                                                              ncols=3,
                                                              row_lim=5,
                                                              figure=figure)
        ax = {}
        ax['drag'] = ax_drag[0]
        ax['rdrag'] = ax_drag[1]
        ax['dis_l'] = ax_dissens[0]
        ax['dis_r'] = ax_dissens[0].twinx()
        ax['rdis'] = ax_dissens[1]
        ax['trap_l'] = ax_trap[0]
        ax['trap_r'] = ax_trap[0].twinx()
        ax['rtrap'] = ax_trap[1]

        self.plot_drag(axis=ax['drag'], height_unit=height_unit)
        self.plot_drag(axis=ax['drag'],
                       plot_outliers=True,
                       height_unit=height_unit)
        self.plot_rel_drag_fit(height_fit=True, axis=ax['drag'])

        self.plot_rel_drag_fit(height_fit=True,
                               plot_residuals=True,
                               axis=ax['rdrag'])

        zplotted = False

        for idx, name in enumerate(names):
            if name.find('z') < 0:
                lor = '_l'
            else:
                # usually the axial values are far off from the lateral values
                # so let's plot them against the right y-axis
                zplotted = True
                lor = '_r'

            # Data + fit
            self.plot_dissens(names=name,
                              axis=ax['dis' + lor],
                              plot_ac_data=False,
                              plot_corrected=plot_corrected,
                              height_unit=height_unit,
                              dissens_unit=dissens_unit,
                              color=col_dict[idx])
            self.plot_dissens(names=name,
                              axis=ax['dis' + lor],
                              plot_ac_data=False,
                              plot_corrected=plot_corrected,
                              plot_outliers=True,
                              height_unit=height_unit,
                              dissens_unit=dissens_unit,
                              color=col_dict[idx + 10])
            self.plot_dissens_fit(names=name,
                                  plot_corrected=plot_corrected,
                                  axis=ax['dis' + lor],
                                  height_unit=height_unit,
                                  dissens_unit=dissens_unit,
                                  color=col_dict[idx])
            # Residuals
            self.plot_dissens_fit(names=name,
                                  plot_residuals=True,
                                  axis=ax['rdis'],
                                  height_unit=height_unit,
                                  dissens_unit=dissens_unit)

            # Data + fit
            self.plot_trap_stiffness(names=name,
                                     axis=ax['trap' + lor],
                                     plot_ac_data=False,
                                     plot_corrected=plot_corrected,
                                     height_unit=height_unit,
                                     trap_stiffness_unit=trap_stiffness_unit,
                                     color=col_dict[idx])
            self.plot_trap_stiffness(names=name,
                                     axis=ax['trap' + lor],
                                     plot_ac_data=False,
                                     plot_corrected=plot_corrected,
                                     plot_outliers=True,
                                     height_unit=height_unit,
                                     trap_stiffness_unit=trap_stiffness_unit,
                                     color=col_dict[idx + 10])
            self.plot_trap_stiffness_fit(names=name,
                                         plot_corrected=plot_corrected,
                                         axis=ax['trap' + lor],
                                         height_unit=height_unit,
                                         trap_stiffness_unit=trap_stiffness_unit,
                                         color=col_dict[idx])

            self.plot_trap_stiffness_fit(names=name,
                                         plot_residuals=True,
                                         axis=ax['rtrap'],
                                         height_unit=height_unit,
                                         trap_stiffness_unit=trap_stiffness_unit,
                                         color=col_dict[idx])
        # do not show the x-axis label for the upper plots
        ax_noxlabel = [a for n, a in ax.items() if not(n.startswith('r'))]
        plt.setp(ax_noxlabel, xlabel='')

        ax_notx = [a.get_xticklabels()
                   for n, a in ax.items()
                   if not(n.startswith('r'))]
        plt.setp(ax_notx, visible=False)
        if not zplotted:
            ax_noylabel = [a for n, a in ax.items()
                           if (not(n.startswith('r')) and n.endswith('r'))]
            plt.setp(ax_noylabel, ylabel='')

        figure.tight_layout()

        if save:
            if path is None:
                flnm = self.basename + '_results_plot.' + save_as
                path = join(self.directory, flnm)
            figure.patch.set_alpha(1.0)
            figure.patch.set_color('white')
            figure.savefig(path, dpi=dpi, format=save_as)
        return figure


class HeightCalibTime(object):
    """
    Use a PyOTI's Motion object to create a HeightCalibrtion object.

    Workflow:
    - __init__(motion, ex_freq, ..)
    - set_exp_conditions(temperature, radius, material='..', medium='..', temp_unit='K', ...)
    - gen_height_calibration(N_avg, pos_signal_unit='um', invert_z_signal=True, ...)
    """
    def __init__(self, motion, ex_freq, psd_names=None, position_names=None):
        """
        Initializes a HeightCalibTime object from a PyOTI Motion object.
        
        Arguments
        ---------
        motion : Motion
        ex_freq : float
            Excitation frequency in Hz.
        psd_names : list of str
            List of names of the data of which to create power spectra from
            (Defaults to ['psdX', 'psdY', 'psdZ']).
        position_names : list of str
            Names of the stage monitor signals.

        Notes
        -----
        names of traces
            Note that the names are either defined in the setup specific
            configuration file or the traces argument when creating a record.
        """
        self.psd_names = psd_names or ['psdX', 'psdY', 'psdZ']
        self.pos_names = position_names or ['positionX', 'positionY', 'positionZ']
        self.motion = motion
        self.expset = None
        self.ex_freq = ex_freq

        for name in self.psd_names:
            if name not in self.traces:
                raise Exception('psd name "{}" not found in traces.'.format(name))
        for name in self.pos_names:
            if name not in self.traces:
                raise Exception('position name "{}" not found in traces.'.format(name))

    @property
    def f_sample(self):
        return self.motion.resolution

    @property
    def ex_psd_axis(self):
        return self.motion.region._excited(traces=['psdX', 'psdY'], index=False)

    @property
    def ex_pos_axis(self):
        return self.motion.region._excited(index=False)

    @property
    def get_data(self):
        return self.motion.get_data

    @property
    def traces(self):
        return self.motion.region.traces

    def set_exp_conditions(self,
                           temp,
                           radius,
                           **kwargs
                           ):
        """
        Set up the experimental conditions.

        Arguments
        ---------
        temp : Temperature in temp_unit (default is kelvin)
        radius : Particle diameter in radius_unit (default is meter)

        Keyword Arguments
        -----------------
        kwargs are passed to pyotc.ExpSetting.__init__().
        defaults :
            temp_err=0
            radius_err=0
            height=inf
            density_particle=None
            density_medium=None
            viscosity=None
            material=''
            medium='water'
            temp_unit='K'
            radius_unit='m'
            height_unit='m'
            density_particle_unit='kg/m**3'
            density_medium_unit='kg/m**3'
            viscosity_unit='Pa*s'
            warn=True
        """
        # define experimental parameters
        self.expset = ExpSetting(temp, radius, **kwargs)

    def gen_height_calibration(self,
                               N_avg,
                               pos_signal_unit='um',
                               invert_z_signal=True,
                               ac_calibration=True,
                               verbose=False,
                               kws_psdfit={}):
        """
        """
        hc = HeightCalibration()
        psd_data = {name: self.get_data(name).flatten()
                    for name in self.psd_names}
        ex_pos_data = self.get_data(self.ex_pos_axis).flatten()
        zpos = self.get_data(self.motion.traces_sf).flatten()

        for plateau in self.motion.plateaus:
            height = (-1 if invert_z_signal else 1) * zpos[plateau].mean()
            es = self.expset.copy()
            es.set_height(height, pos_signal_unit)

            pm = PSDMeasurement(exp_setting=es)

            for psd_name, psd_signal in psd_data.items():
                if verbose:
                    print('generating psd {0} at {1:1.3f} {2}'
                          ''.format(psd_name, height, pos_signal_unit))

                p = gen_PSD_from_time_series(psd_signal[plateau],
                                             self.f_sample,
                                             N_avg,
                                             name=psd_name)

                pm.add_psd(psd_name, p)

                # determine the power in the peak
                if ac_calibration and psd_name is self.ex_psd_axis:
                    if verbose:
                        print('determining power in peak of psd {}'
                              ''.format(psd_name))
                    pp = gen_PSD_from_time_series(psd_signal[plateau],
                                                  self.f_sample,
                                                  N_avg,
                                                  calc_errors=True)
                    ex_power = pp.psd[pp.freq == self.ex_freq] * pp.df
                    ex_power_err = pp.psd_err[pp.freq == self.ex_freq] * pp.df

            if ac_calibration:
                # get exciation amplitude from stage signal
                if verbose:
                    print('determining power in stage signal')
                p = gen_PSD_from_time_series(ex_pos_data[plateau],
                                             self.f_sample,
                                             N_avg,
                                             calc_errors=True)
                # amplitude = sqrt(2) * A_rms = sqrt(2) * sqrt(P(f_drive) * df)
                ex_amplitude = sqrt(2 * p.psd[p.freq == self.ex_freq] * p.df)
                ex_amplitude_err = sqrt(2 * p.psd_err[p.freq == self.ex_freq]
                                        * p.df)

                pm.set_ac_params(self.ex_psd_axis,
                                 self.ex_freq,
                                 ex_amplitude[0],
                                 ex_power[0],
                                 ex_amplitude_err=ex_amplitude_err[0],
                                 ex_power_err=ex_power_err[0],
                                 amplitude_unit=pos_signal_unit)

            pf = PSDFit(pm, **kws_psdfit)
            hc.add_psdfit(height, pf, unit=pos_signal_unit)
            if verbose:
                print('created fit object and added to HC-object.')
            else:
                print('.', end='')

        return hc
