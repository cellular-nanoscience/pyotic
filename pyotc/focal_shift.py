# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
"""
Provides functions to find the focal shift of an optical tweezers setup.
"""

from . import ureg

from lmfit import minimize
from lmfit import Parameters
from lmfit import report_fit

import matplotlib.pyplot as plt

from scipy import exp
from scipy import mean
from scipy import pi
from scipy import sin

from .plotting import add_plot_to_figure


def fs_fun(params, heights, wavelength, ref_index, signal=None, simple=False):
    """
    Calculates the difference (residual) between the predicted intensity
    from the model function and the experimental data.

    Arguments
    ---------
    params : lmfit.Parameters
        Object that hold "lmfit.Parameter()"
        Parameters are:
            focal_shift, U, tau, h0, p0, .., p3
    height : array(float)
        Apparent distance between the surface and the bead center in um.
    wavelength : float
        The wavelength of the trapping laser in nm.
    ref_index : float
        Refractive index of the sample medium.
    signal : array(float) or None
        Intensity-signal to calculate the residuals to. If None, the model for
        the given parameters is returned.
    """
    h = heights
    p = params.valuesdict()
    fs = p['focal_shift']
    U = p['U']
    tau = p['tau']
    h0 = p['h0']
    p0 = p['p0']
    if not simple:
        p1 = p['p1']
        p2 = p['p2']
        p3 = p['p3']

    wl = wavelength * 1e-3 / (2 * ref_index)

    if simple:
        model = (U * exp(-(h - h0)**2 / tau) *
                 sin((2 * pi / wl) * fs * (h - h0)) +
                 p0)
    else:
        model = (U * exp(-(h - h0)**2 / tau) *
                 sin((2 * pi / wl) * fs * (h - h0)) +
                 p0 + p1 * h + p2 * h**2 + p3 * h**3)

    if signal is None:
        return model

    return (signal - model)


def gen_fs_fit_pars(focal_shift=0.8,
                    U=0.5,  # Volts
                    h0=0.0,  # stage offset in um
                    tau=1.0,  # decay rate of gaussian slope
                    p0=8.5,  # mean intensity signal
                    p1=0.01,  # polynomial parameters
                    p2=0.01,
                    p3=0.001,
                    simple=False
                    ):
    """ Generate the lmfit fitting Parameters object. """
    pars = Parameters()
    pars.add('focal_shift', value=focal_shift, min=0, max=1.0)
    pars.add('U', value=U)
    pars.add('h0', value=h0)
    pars.add('tau', value=tau, min=0)
    pars.add('p0', value=p0)
    if not simple:
        pars.add('p1', value=p1)
        pars.add('p2', value=p2)
        pars.add('p3', value=p3)

    return pars


def get_focal_shift(heights,  # in um
                    signal,
                    wavelength,  # in nm
                    ref_index,
                    simple=False,
                    height_unit='um',
                    wavelength_unit='nm',
                    parameters=None,
                    report_focal_shift=False,
                    plot_fit=False,
                    plt_axis=None,
                    verbose=False,
                    **init_kwargs
                    ):
    """
    Find the focal shift of the optical trap that is caused due to the
    refravtive index mismatch between coverslip and medium.

    The function uses teh following model to find the focal shift

    Model
        I = (U * exp(-(h - h0)**2 / tau) *
             sin(2 * pi * (h - h0) / wl) +
             p0 + p1 * h + p2 * h**2 + p3 * h**3)
        where wl = wavelength * 1e-3 / (2 * ref_index * fs)

    Arguments
    ---------
    heights : array
        Distances between surface and bead-center.
    signal : array
        Intensities that show interference fluctuations, due to standing
        wave between surface and trapped bead.
    wavelength : float
        Wavelength of the trapping laser.
    ref_index : float
        Refractive index of the medium.
    simple : bool
        If True, the function only uses an offset p0, if False, the function
        uses a 3rd-order polynomial.
    height_unit : str
        Unit of the heights.
    wavelength_unit : str
        Unit of the wavelength.
    parameters : lmfit.Parameters()
        Object with predefined initial guesses of the fit function.
        If None predefined values are used. Use gen_fs_fit_pars() to
        initialize.
    verbose : bool
        Plot fit and give more info if True.

    Keyword Arguments
    -----------------
    focal_shift
    U
    h0
    tau
    p0
    p1
    p2
    p3

    Returns
    -------
    minimizer : lmfit.Minimizer

    See Also
    --------
    gen_fs_fit_pars
    """
    if height_unit != 'um':
        heights *= ureg(height_unit).to('um').magnitude
    if wavelength_unit != 'nm':
        wavelength *= ureg(wavelength_unit).to('nm').magnitude

    # define parameters for the fit
    if parameters is None:
        init_pars = {'focal_shift': 0.8,
                     'U': max(signal)-min(signal),
                     'h0': min(heights)+0.2*(min(heights)+max(heights)),
                     'p0': mean(signal)
                     }
        for k, v in init_pars.items():
            if k not in init_kwargs:
                init_kwargs[k] = v
        pars = gen_fs_fit_pars(simple=simple, **init_kwargs)
    else:
        pars = parameters

    minimizer = minimize(fs_fun,
                         pars,
                         args=(heights,  # in um
                               wavelength,  # in nm
                               ref_index),
                         kws={'signal': signal,
                              'simple': simple})

    minimizer.fcn_eval = fs_fun(minimizer.params,
                                heights,
                                wavelength,
                                ref_index,
                                simple=simple)
    if verbose:
        report_fit(minimizer)
    if plot_fit or verbose:
        if plt_axis:
            fig = plt_axis.figure
        else:
            fig = plt.figure()
        ax = add_plot_to_figure(fig,
                                heights, signal,
                                fmt='o',
                                axis=plt_axis)

        add_plot_to_figure(fig, heights, minimizer.fcn_eval, fmt='-', axis=ax)
        fig.show()
    if minimizer.success:
        fs = minimizer.params['focal_shift'].value
        dfs = minimizer.params['focal_shift'].stderr
        if report_focal_shift or verbose:
            print('Focal shift = {0:1.3f} +/- {1:1.3f}'.format(fs, dfs))
    else:
        if verbose:
            print('Least squares fit did not succeed.')
            print('Minimizer message: "{0:s}"'.format(minimizer.message))
            print('lmdif message: "{0:s}"'.format(minimizer.lmdif_message))

    return minimizer
