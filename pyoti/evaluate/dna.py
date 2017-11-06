# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:29:39 2016

@author: Tobias Jachowski
"""
import numpy as np
from scipy import constants
from collections import namedtuple

k_B = constants.value('Boltzmann constant')
ForceExtension = namedtuple('ForceExtension', ['extension', 'force'])


def worm_like_chain(x, L_0, L_p=50e-9, T=298):
    """
    Marko, J.F.; Eric D. Siggia (1995). "Stretching DNA". Macromolecules. 28:
    8759–8770. doi:10.1021/ma00130a008

    http://biocurious.com/2006/07/04/wormlike-chains

    Parameters
    ----------
    x : float
        extension (m)
    L_0 : float
        contour length (m)
    L_p : float
        persistence length (m)
    T : float
        temperature (K)
    """
    F = k_B * T / L_p * (1 / (4 * (1 - x / L_0)**2) - 1/4 + x / L_0)

    return F


def worm_like_chain_1p(x, L_0, L_p=50e-9, T=298):
    """
    Petrosyan, R. (2016). "Improved approximations for some polymer extension
    models". Rehol Acta. doi:10.1007/s00397-016-0977-9

    Parameters
    ----------
    x : float
        extension (m)
    L_0 : float
        contour length (m)
    L_p : float
        persistence length (m)
    T : float
        temperature (K)
    """
    F = k_B * T / L_p * (1 / (4 * (1 - x / L_0)**2) - 1/4 + x / L_0
                         - 0.8 * (x / L_0)**2.15)

    return F


# TODO: Check units
def twistable_wlc(F, L_0, x=None, L_p=43.3e-9, K_0=1246e-12, S=1500e-12,
                  C=440e-30, T=298.2):
    """
    Twistable worm like chain model Gross et al. 2011

    Paramaters
    ----------
    F : float
        force (N)
    L_0 : float
        contour length (m)
    x : float
        extension (m)
    L_p : float
        persistence length (m)
    K_0 : float
        elastic modulus (N)
    S : float
        stretch modulus (N)
    C : float
        twist rigidity (Nm²)
    """
    def g(F, fitting=False, use_lambda_fit=False, use_unwinding_fit=False):
        """
        g(F) describes the twist-stretch coupling how DNA complies to tension.

        The original equation is as follows:
        g(F) = (S * C - C * F * (x/L_0 - 1 + 1/2
                * (k*T/(F*L_p))**(1/2))**(-1))**(1/2)

        The equation can be simplified:
            - below Fc of 30 pN: as a constant (- 100 pN nm)
            - above Fc to linear order: g0 + g1 * F

        g0 was determied to be 590 pN nm (or 560 pN nm with lambda DNA)
        g1 was determined to be 18 nm
        """
        if fitting:
            return (S * C - C * F * (x/L_0 - 1 + 1/2
                    * (k_B*T/(F*L_p))**(1/2))**(-1))**(1/2)
        if F <= 30e-12:  # N
            return - 100e-21  # Nm
        else:
            g0 = - 590e-21  # Nm
            if use_lambda_fit:
                g0 = - 560e-21  # Nm
            if use_unwinding_fit:
                g0 = - 637e-21  # Nm
            return g0 + 17e-9 * F

    x = L_0 * (1 - 1/2 * (k_B * T / (F * L_p))**(1/2)
               + C / (-g(F)**2 + S * C) * F)

    return x


def force_extension(bps, pitch=0.338e-9, L_p=43.3e-9, T=298.2, min_ext=0.5,
                    max_ext=0.978, samples=1000):
    """
    Parameters
    ----------
    bps : int
        number of base-pairs
    pitch : float
        length of one base-pair (m)
        Defaults to 0.338e-9 m (Saenger 1988)
    L_p : float
        persistence length (m)
    T : float
        temperature (K)
    min_ext : float
        minimum extension, normalized to contour length
        0.5-0.978 entspricht ~520nm-1057nm (0.118...49.2pN)
    max_ext : float
        maximum extension, normalized to contour length
    samples : int
        number of samples to generate

    Returns
    -------
    ForceExtension
        The namedtuple consists of extension and force in (m)
    """
    L_0 = bps * pitch  # m contour length
    x = np.linspace(min_ext * L_0, max_ext * L_0, samples)
    F = worm_like_chain(x, L_0, L_p, T)

    # K_0=1246e-12  # N elastic modulus
    # x = (x + f / K_0) * L_0

    return ForceExtension(extension=x, force=F)


def twistable_force_extension(bps, pitch=0.338e-9, L_p=43.3e-9, K_0=1246e-12,
                              S=1500e-12, C=440e-30, T=298.2, min_f=10e-12,
                              max_f=60e-12, samples=1000):
    F = np.linspace(min_f, max_f, samples)
    L_0 = bps * pitch
    x = twistable_wlc(F, L_0, L_p, K_0, S, C, T)

    return ForceExtension(extension=x, force=F)


import operator
def crop_x_y(x, y=None, min_x=None, max_x=None, min_y=None, max_y=None,
             include_bounds=True):
    """
    Crop pairs of variates according to their minimum and maximum values.

    Parameters
    ----------
    x : 1D numpy.ndarray of type float
        The x values.
    y : 1D numpy.ndarray of type float
        The y values.
    min_x : float
        The minimum value of `x`.
    max_x : float
        The maximum value of `x`.
    min_y : float
        The minimum value of `y`.
    max_y : float
        The maximum value of `y`.

    Returns
    -------
    tuple of 2 1D numpy.ndarray of type float
        The cropped values (x, y).
    """
    if include_bounds:
        ol = operator.le
        og = operator.ge
    else:
        ol = operator.lt
        og = operator.gt
    max_x = max_x or float('inf')
    min_x = min_x or float('-inf')
    i_x = ol(x, max_x)
    i_x = np.logical_and(i_x, og(x, min_x))
    if y is None:
        return x[i_x]
    max_y = max_y or float('inf')
    min_y = min_y or float('-inf')
    i_y = ol(y, max_y)
    i_y = np.logical_and(i_y, og(y, min_y))
    i = np.logical_and(i_x, i_y)
    return (x[i], y[i])


def residual(params, model_func, x, data, eps=None):
    """
    Calculate the residuals of a model and given data.

    Parameters
    ----------
    fit_func : function
    x : 1D numpy.ndarray of type float
    params : dict
    data : 1D numpy.ndarray of type float
    eps : float
    """
    model = model_func(x, **params)

    # Calculate the residuals of the model and the given data
    if eps is None:
        return model - data
    return (model - data) / eps


def crop_fe_wlc_dna(e, f, L_0, min_e=None, max_e=None, max_f=None):
    # Choose boundaries to avoid nan values and max force of 25 pN, up to where
    # wlc is valid
    min_x = min_e or 0.01e-9
    max_x = max_e or L_0
    max_f = max_f or 15e-12
    return crop_x_y(e, f, min_x=min_x, max_x=max_x, max_y=max_f,
                    include_bounds=False)


def residual_wlc_dna(params, model_func, e, f, min_e=None, max_e=None,
                     max_f=None):
    # Crop data to be fitted
    #   a) avoid nan values, and
    #   b) limit to force where model is valid
    L_0 = params['L_0']
    _e, _f = crop_fe_wlc_dna(e, f, L_0, min_e, max_e, max_f)
    return residual(params, model_func, _e, _f)


from lmfit import minimize, Parameters, fit_report
def fit_force_extension(e, f, bps, pitch=0.338e-9, L_p=50e-9, T=298,
                        model_func=None, fix_L_0=False, fix_L_p=False,
                        fix_T=True, min_e=None, max_e=None, max_f=None,
                        verbose=False):
    # Calculate the contour length of the DNA
    L_0 = bps * pitch  # m contour length

    params = Parameters()
    params.add('L_0', value=L_0, vary=not fix_L_0)
    params.add('L_p', value=L_p, vary=not fix_L_p)
    params.add('T', value=T, vary=not fix_T)

    model_func = model_func or worm_like_chain

    # Do the fitting ...
    out = minimize(residual_wlc_dna, params, args=(model_func, e, f, min_e,
                                                   max_e, max_f))

    if verbose:
        print(fit_report(out))
        print('[[DNA related info]]')
        print('    Number of base-pairs: {:.0f}'.format(
                                        np.round(out.params['L_0'] / pitch)))

    return out
