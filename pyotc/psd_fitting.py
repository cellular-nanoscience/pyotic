# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
from . import ureg

from . import name_constants as co

from .plotting import add_plot_to_figure
from .plotting import col_dict

from .psd import apply_aliasing
from .psd import apply_low_pass_filter
from .psd import lorentzian_psd
from .psd import hydro_psd

from .utilities import check_and_separate
from .utilities import flatten_list
from .utilities import gen_fit_pars
from .utilities import logspace_points_per_decade

from collections import OrderedDict

from configparser import ConfigParser

from functools import partial

from lmfit import minimize
from lmfit import Model
from lmfit import report_fit

import matplotlib.pyplot as plt

from os.path import join, isfile

from scipy import array
from scipy import complex_
from scipy import inf
from scipy import isnan
from scipy import logical_and
from scipy import pi
from scipy import sqrt

from scipy.constants import Boltzmann as k_B
from scipy.stats import norm

import sys

import time

import warnings


def gen_psd_fit_pars(D=1, D_min=0,
                     f_c=1000, f_c_min=0, f_c_max=1e7,
                     **kwargs):
    """
    Generate psd fit parameters.

    Every keyword will create a parameter, keywords ending with '_min' or
    '_max' will define the upper and lower bound for the parameter.

    Example
    -------
    gen_fit_pars(D=1e-6, D_min=0, f_c=1000, f_c_max=1e6) will create
    OrderedDict([('D', <Parameter 'D', 1e-06, bounds=[0:None]>),
    ('f_c', <Parameter 'f_c', 1000, bounds=[None:1000000.0]>)]).

    Returns
    -------
    params : Parameters
        lmfit Parameters object
    """
    params = gen_fit_pars(D=D, D_min=D_min,
                          f_c=f_c, f_c_min=f_c_min, f_c_max=f_c_max,
                          **kwargs)
    return params


def gen_model_fun(model='lorentzian',
                  lp_filter=False,
                  lp_fixed=False,
                  aliasing=False,
                  f_sample=None,
                  N_alias=9,
                  model_fun=None,
                  param_names=None,
                  name='other',
                  invert_function=False,
                  **kwargs):
    """
    Generate the model function.

    invert function has no effect on userdefined model_functions.
    They must be inverted already!
    """
    ################################
    # get the basic model function #
    ################################
    model = model.lower()
    if model == 'lorentzian':
        model_fun = lorentzian_psd
        name = 'Lorentzian PSD'
    elif model == 'hydro':
        model_fun = hydro_psd
        name = 'Hydr. PSD'
    elif model == 'other':
        if model_fun is None:
            raise Exception('Function missing in keyword argument model_fun')
    else:
        raise Exception('Model unknown model={}'.format(model))

    ###################################
    # define the fitting parameters #
    ###################################
    if param_names is None:
        param_names = ['D', 'f_c']
        if lp_filter:
            name += ', low-pass filtered'
            param_names.append('f3dB')
            param_names.append('alpha')

    # get rid of init guess kwargs, which are not important here
    for pname in param_names:
        try:
            kwargs.pop(pname)
        except:
            pass
    #################################
    # generate the fitting function #
    #################################
    # separate the function keyword arguments from all keyword arguments
    model_kwargs = check_and_separate(model_fun, kwargs, exclude=True)

    # basic model function (partial), e.g. hydro_psd(freq, D, f_c)
    pfun = partial(model_fun, **model_kwargs)
    pfun.__name__ = model_fun.__name__
    if not model_fun.__doc__:
        model_fun.__doc__ = ''
    pfun.__doc__ = ('partial of {0:s}\nOriginal doc\n{1:s}'
                    ''.format(pfun.__name__, model_fun.__doc__))

    # if user defined function is provided don't do aliasong etc.
    # this must have been don by the user beforehand
    if model != 'other':
        # define function w/ resp. to lp_filter and aliasing and INVERT
        if aliasing:
            name += ' and aliased'

            if f_sample is None:
                raise Exception('Sampling frequency is missing.')

            if lp_filter:
                def aliased_lp_fun(freq, D, f_c, f3dB, alpha):
                    f = apply_low_pass_filter(pfun, f3dB, alpha)
                    f_ = apply_aliasing(f, f_sample, N_alias)(freq, D, f_c)
                    if invert_function:
                        return 1 / f_
                    else:
                        return f_

                model_fun_ = aliased_lp_fun
            else:
                def aliased_fun(freq, D, f_c):
                    f = apply_aliasing(pfun, f_sample, N_alias)(freq, D, f_c)
                    if invert_function:
                        return 1 / f
                    else:
                        return f

                model_fun_ = aliased_fun
        else:
            if lp_filter:
                def lp_fun(freq, D, f_c, f3dB, alpha):
                    f = apply_low_pass_filter(pfun, f3dB, alpha)(freq, D, f_c)
                    if invert_function:
                        return 1 / f
                    else:
                        return f

                model_fun_ = lp_fun
            else:
                def fun(freq, D,  f_c):
                    if invert_function:
                        return 1 / pfun(freq, D, f_c)
                    else:
                        return pfun(freq, D, f_c)

                model_fun_ = fun
    else:
        model_fun_ = pfun

    return (model_fun_, param_names, name)


def fit_psd(freq,
            psd_data,
            psd_err=None,
            N_avg=1,
            model='lorentzian',
            lp_filter=False,
            lp_fixed=False,
            aliasing=False,
            f_sample=None,
            N_alias=9,
            model_fun=None,
            name='other',
            verbose=False,
            **kwargs):
    """
    Fits the given power spectral density to the selected model.

    kwargs holds both model_kwargs and inital guess and min max kwargs.
    """
    #################################
    # generate the fitting function #
    #################################
    # kwargs still holds both, model_kwargs and inital guesses, min max etc.
    model_fun_, param_names, name = gen_model_fun(model=model,
                                                  invert_function=True,
                                                  lp_filter=lp_filter,
                                                  lp_fixed=lp_fixed,
                                                  aliasing=aliasing,
                                                  f_sample=f_sample,
                                                  N_alias=N_alias,
                                                  model_fun=model_fun,
                                                  name=name,
                                                  **kwargs)
    # kwargs should only have the initial guesses and min, max etc.

    ###################################
    # generate the fitting parameters #
    ###################################
    # get initial guesses and mins and maxs from kwargs or create them
    par_kwargs = {}
    suffix = ['', '_min', '_max']
    for pname in param_names:
        for suff in suffix:
            try:
                par_kwargs[pname + suff] = kwargs.pop(pname + suff)
            except:
                pass
    if lp_filter:
        if 'f3dB' not in par_kwargs:
            par_kwargs['f3dB'] = 10e3
            if lp_fixed:
                warnings.warn('cut-off frequency not provided! '
                              'Using 10 kHz as fallback.')
        if 'f3dB_min' not in par_kwargs:
            par_kwargs['f3dB_min'] = 0

        if 'alpha' not in par_kwargs:
            par_kwargs['alpha'] = 0.5
            if lp_fixed:
                warnings.warn('alpha not provided! Using 0.5 as fallback.')
        if 'alpha_min' not in par_kwargs:
            par_kwargs['alpha_min'] = 0
        if 'alpha_max' not in par_kwargs:
            par_kwargs['alpha_max'] = 1

    # generate the Parameters object
    pars = gen_psd_fit_pars(**par_kwargs)

    if lp_fixed:
        pars['f3dB'].vary = False
        pars['alpha'].vary = False

    ############################
    # generate the lmfit.Model #
    ############################
    model_obj = Model(model_fun_,
                      independent_vars=['freq'],
                      param_names=param_names,
                      name=name)

    #################################
    # calculate the weights and fit #
    #################################
    if psd_err is None:
        # err = psd_data / sqrt(N_avg)
        # err_i = 1 / (psd_data * sqrt(N_avg))
        weights = psd_data * sqrt(N_avg)
    else:
        # err = psd_data / sqrt(N_avg)
        # err_i = 1 / (psd_data * sqrt(N_avg)) = 1 / (psd_err * N_avg)
        # err_i = 1 / (err * N_avg)
        weights = psd_err * N_avg

    modelfit = model_obj.fit(1/psd_data, freq=freq,
                             params=pars, weights=weights,
                             fit_kws={'nan_policy': 'propagate'})
    if verbose:
        report_fit(modelfit)

    return modelfit


def analytical_lsq_lorentzian(freq, psd, verbose=False):
    """
    Analytically solve the least-squares problem of a Lorentzian-shaped power
    spectrum.

    Calculates the least-squares fit of a Lorentzian power spectrum
    analytically for the given frequencies and power spectral density values
    (see Equs. 13 and 14 in Berg-Sørensen & Flyvbjerg (2004))

    Arguments
    ---------
    freq : array(float)
        Frequency vector.
    psd : array(float)
        psd values.
    verbose : bool
        Will give additional information, if True.

    Returns
    -------
    f_c : float
        The corner frequency
    D : float
        The diffusion coefficient
    chi2 : float
        The reduced chi-squared value of the fit.

    Note
    ----
    Bias
        A biased diffusion constant is NOT accounted for as it is stated in
        Equ. (41) in Nørrelykke & Flyvbjerg (2010).

    Eaution
        A factor two is omitted in the formula for D because we are dealing with
        one-sided power spectra.

    References
    ----------
    Berg-Sørensen, K., & Flyvbjerg, H. (2004)
        Equ. (13, 14) in
        Power spectrum analysis for optical tweezers. Review of Scientific
        Instruments, 75(3), 594–612.
        http://doi.org/10.1063/1.1645654

    Nørrelykke & Flyvbjerg (2010)
        Power spectrum analysis with least-squares fitting:
        Amplitude bias and its elimination, with application to optical
        tweezers and atomic force microscope cantilevers.
        Review of Scientific Instruments, 81(7), 075103.
        http://doi.org/10.1063/1.3455217
    """
    def ls_sum(freq, psd, p, q):
        return sum(freq**(2 * p) * psd**q)

    S_00 = ls_sum(freq, psd, 0, 0)
    S_01 = ls_sum(freq, psd, 0, 1)
    S_02 = ls_sum(freq, psd, 0, 2)
    # S_10 = ls_sum(freq, psd, 1, 0)
    S_11 = ls_sum(freq, psd, 1, 1)
    S_12 = ls_sum(freq, psd, 1, 2)
    # S_21 = ls_sum(freq, psd, 2, 1)
    S_22 = ls_sum(freq, psd, 2, 2)

    a = (S_01 * S_22 - S_11 * S_12)
    b = (S_11 * S_02 - S_01 * S_12)

    f_c = sqrt(a / b).real
    D = pi * pi * (S_02 * S_22 - S_12 * S_12) / (S_11 * S_02 - S_01 * S_12)
    chi2 = S_00 - ((S_01 * S_01 * S_22 +
                    S_11 * S_11 * S_02 - 2 * S_01 * S_11 * S_12) /
                   (S_02 * S_22 - S_12 * S_12))

    if isnan(f_c):
        if verbose:
            warnings.warn('Corner frequency calculation resulted in NAN.'
                          'f_c = (({0:1.17e} / {1:1.17e})**0.5).real'
                          ''.format(a, b))
            print('The following sums were calculated:'
                  'S_00={}\t S_01={}\t S_02={}\n'
                  'S_11={}\t S_12={}\t S_22={}'
                  ''.format(S_00, S_01, S_02, S_11, S_12, S_22))
            print('with fmin = {}, fmax={}'.format(min(freq), max(freq)))
            print('a=(S_01 * S_22 - S_11 * S_12)')
            print('b=(S_11 * S_02 - S_01 * S_12)')

    return (f_c, D, chi2)


def calc_anal_lsq_opt(freq, psd, verbose=False, plot=False, ppd=20):
    """
    Try to find the optimal range for the analytical least-squares algorithm,
    by varying the range of the frequency vector.

    The function varies the range of data points to find an optimal least-
    squares fit of the data. this is done as follows:
    First, the highest frequency boundary is varied, while the lower is kept
    fixed. For each new set of lower and upper bound the
    *analytical_lsq_lorentzian*-method is called. The frequency, at which the
    chi2 value is minimal, is takes as upper bound.
    Second, the upper frequency bound is varied and the frequency, at which the
    corner frequency is highest is taken as lower bound.

    Arguments
    ---------
    freq : array(float)
        Frequency vector.
    psd : array(float)
        psd values.
    verbose : bool
        Will give additional information, if True.
    plot : bool
        Whether, to plot the trials (debugging)
    ppd : int
        Number of point per decade. The more the higer the number of trials
        (default is 20).

    Returns
    -------
    f_c : float
        The corner frequency of the *best* fit.
    D : float
        The diffusion coefficient of the *best* fit.
    chi2 : float
        The reduced chi-squared value of the *best* fit.
    (fmin, fmax) : tuple
        Tuple defining the boundaries of the *best* fit.

    See Also
    --------
    analytical_lsq_lorentzian
    """
    if verbose:
        print('Trying to find optimal boundaries for the analytical least '
              'squares of a Lorentzian.')
        print('Number of points per decade: {}'.format(ppd))

    fmin_init = min(freq)

    # decimate frequency points
    #f = logspace(log10(min(freq)),
    #             log10(max(freq)),
    #             base=10, num=n_trials)
    f = logspace_points_per_decade(min(freq), max(freq), ppd=ppd)
    # define a set of test values for fmin and fmax
    fmin_ = f[:int(len(f)/2)]
    fmax_ = f[int(len(f)/2):]

    if verbose:
        print('initial fmin: {0:6.2f}'.format(fmin_init))
        print('fmin - fmax border: {0:6.2f}'.format(f[round(len(f) / 2)]))

    # vary the upper freq bound and get corner frequency and chi² value
    # from the analytical least squares calculation
    fc2 = []
    ch2 = []
    for fmax in fmax_:
        excl = ((freq < fmin_init) + (freq > fmax))
        keep = ~excl

        freq_part = freq[keep]
        res = analytical_lsq_lorentzian(freq_part, psd[keep], verbose=verbose)
        if not(any(isnan(a) for a in res)):
            fc2.append(res[0])
            ch2.append(res[2] / len(freq_part))

    # set the upper bound where chi² is minimal
    if len(ch2) > 0:
        fmax = array(fmax_)[ch2 == min(ch2)][0]
    else:
        fmax = max(freq)
    # plot the corner frequencies and chi² values vs the upper
    # frequency bound
    if plot:
        fig = plt.figure()
        ax1 = add_plot_to_figure(fig, [fmax, fmax], [max(ch2), min(ch2)],
                                 subplot=(1, 3, 1))
        ax1 = add_plot_to_figure(fig, fmax_, fc2,
                                 title=(r'$\mathsf{f_c and \chi^2 vs. '
                                        'f_{max}}$'),
                                 ylabel=r'$\mathsf{f_c (Hz)}$',
                                 fmt='or',
                                 axis=ax1,
                                 label=r'$\mathsf{f_c}$',
                                 showLegend=True,
                                 legend_kwargs={'loc': 2})
        ax1a = ax1.twinx()
        add_plot_to_figure(fig, fmax_, ch2,
                           xlabel=r'$\mathsf{f_{max}}$',
                           ylabel=r'$\mathsf{\chi^2}$',
                           fmt='og', axis=ax1a,
                           label=r'$\mathsf{\chi^2}$',
                           showLegend=True)
        #while not plt.waitforbuttonpress():
        #    pass
        #ax.figure.clf()

    if verbose:
        print('optimal fmax: {0:1.3e}'.format(fmax))

    # now vary the lower frequency bound with fixed ('optimal')
    # upper  bound. Calculate again the analytical LSq
    fc1 = []
    ch1 = []
    for fmin in fmin_:
        excl = ((freq < fmin) + (freq > fmax))
        keep = ~excl

        freq_part = freq[keep]
        res = analytical_lsq_lorentzian(freq_part, psd[keep])

        if not(any(isnan(a) for a in res)):
            fc1.append(res[0])
            ch1.append(res[2] / len(freq_part))
        else:
            fc1.append(-1.0)
            ch1.append(-1.0)

    # take the lower bound where the corner freq is highest
    if len(fc1) > 0:
        fmin = array(fmin_)[fc1 == max(fc1)][0]
    else:
        fmin = fmin_init
    # plot the corner frequencies and the chi² value vs lower freq bound
    if plot:
        fig = plt.figure()
        ax2 = add_plot_to_figure(fig, [fmin, fmin], [max(fc1), min(fc1)],
                                 subplot=(1, 3, 2))
        add_plot_to_figure(fig, fmin_, fc1,
                           ylabel=r'$\mathsf{f_c (Hz)}$',
                           fmt='or',
                           axis=ax2,
                           label=r'$\mathsf{f_c}$',
                           legend_kwargs={'loc': 2},
                           showLegend=True)
        ax2a = ax2.twinx()
        add_plot_to_figure(fig, fmin_, ch1,
                           title=r'$\mathsf{f_c and \chi^2 vs. f_{min}}$',
                           xlabel=r'$\mathsf{f_{min}}$',
                           ylabel=r'$\mathsf{\chi^2}$',
                           fmt='og', axis=ax2a,
                           label=r'$\mathsf{\chi^2}$',
                           showLegend=True)
        #while not plt.waitforbuttonpress():
        #    pass
        #ax.figure.clf()

    if verbose:
        print('optimal fmin: {0:6.2f}'.format(fmin))

    excl = ((freq < fmin) + (freq > fmax))
    keep = ~excl
    freq_part = freq[keep]
    f_c, D, chi2 = analytical_lsq_lorentzian(freq_part, psd[keep])

    if verbose:
        print('(f_c, D, chi2) = ({0:1.2f}, {1:1.3e}, {2:1.2f})'
              ''.format(f_c, D, chi2))

    if plot:
        fig = plt.figure()
        ax3 = add_plot_to_figure(fig, freq, psd,
                                 label='PSD',
                                 fmt='+',
                                 alpha=0.3,
                                 subplot=(1, 3, 3)
                                 )
        analpsd = lorentzian_psd(freq_part, D, f_c)
        add_plot_to_figure(fig, freq_part, analpsd,
                           axis=ax3,
                           title='PSD and analytical Lorentzian fit',
                           xlabel=r'Frequency',
                           ylabel=r'$\mathsf{PSD \; (V^2/Hz)}$',
                           label='Lorentzian fit',
                           fmt='--r',
                           showLegend=True,
                           logplot=True
                           )
    return (f_c, D, chi2, (fmin, fmax))


def collective_psd_fit_fun(params,
                           names,
                           freq_dict,
                           model_fun,
                           param_names,
                           psd_dict=None,
                           err_dict=None,
                           verbose=False
                           ):
    """
    residual function for teh collective psd fit.
    """
    psd_ = []
    err_ = []
    model_ = []

    for name in names:
        kws = {pname: params[pname + '_' + name].value
               for pname in param_names if pname not in ['f3dB', 'alpha']}
        kws['f3dB'] = params['f3dB'].value
        kws['alpha'] = params['alpha'].value

        model_.append(list(model_fun(freq_dict[name], **kws)))

        if psd_dict is not None:
            psd_.append(list(psd_dict[name]))
        if err_dict is not None:
            err_.append(list(err_dict[name]))

    model_flat_array = array(list(flatten_list(model_)))

    if psd_dict is not None:
        psd_flat_array = array(list(flatten_list(psd_)))

    if err_dict is not None:
        err_flat_array = array(list(flatten_list(err_)))

    if psd_dict is None:
        return model_flat_array
    elif err_dict is None:
        return (psd_flat_array - model_flat_array)
    else:
        return ((psd_flat_array - model_flat_array) / err_flat_array)


def make_collective_psd_fit(freq_dict,
                            psd_dict,
                            err_dict,
                            model='lorentzian',
                            aliasing=False,
                            f_sample=None,
                            N_alias=9,
                            verbose=False,
                            **kwargs):
    """
    Fit the provided PSDs collectively, providing the ability to fit low-pass
    filter parameters (alpha and f3dB ) same for all PSDs.
    Returns <lmfit.Minimizer>

    Arguments
    ---------
    freq_dict : dict
        Dictionary of frequency vectors for psds with name "name".
    psd_dict : dict
        Dictionary with names as keys and psds as values.
        If None, the function returns the function values instead of the
        residuals.
    err_dict : dict
        Dictionary with names as keys and psds as values.
        If None, the function returns unweighted residuals, i.e. physical
        errors.
    hydro_kwargs : dict
        Dictionary with keyword-arguments passed to hydro_psd(). I.e.
        {'R': radius, 'height': height, 'temp': Temperature,
        'rho': mass_density}
    params : lmfit.Parameters
        Object of the lmfit.Parameters class, holding Paramerter-objects with
        names 'D_name' and 'f_c_name' and the parameters for the low-pass
        filter 'f3dB' and 'alpha'. If None, they are created.
    init_vals : dict
        Dictionary with 'names' as keys and tuples (D, f_c) of initial
        guess-values for the respective psd fit.
        If None D = 1 and f_c = 1000 Hz is used.

    Return
    ------
    lmfit.Minimizer
    """
    names = psd_dict.keys()

    model_fun, param_names, name = gen_model_fun(model=model,
                                                 lp_filter=True,
                                                 lp_fixed=False,
                                                 aliasing=aliasing,
                                                 f_sample=f_sample,
                                                 N_alias=N_alias,
                                                 **kwargs)

    # get initial guesses and mins and maxs from kwargs or create them
    par_kwargs = {}
    suffix = ['', '_min', '_max']
    for pname in param_names:
        if pname in ['f3dB', 'alpha']:
            for suff in suffix:
                try:
                    par_kwargs[pname + suff] = kwargs.pop(pname + name)
                except:
                    pass
        else:
            for name in names:
                for suff in suffix:
                    try:
                        key = kwargs.pop(pname + '_' + name)
                        par_kwargs[pname + '_' + name + suff] = key
                    except:
                        if suff == '':
                            par_kwargs[pname + '_' + name] = 1

    if 'f3dB' not in par_kwargs:
        par_kwargs['f3dB'] = 10e3
    if 'f3dB_min' not in par_kwargs:
        par_kwargs['f3dB_min'] = 0

    if 'alpha' not in par_kwargs:
        par_kwargs['alpha'] = 0.5
    if 'alpha_min' not in par_kwargs:
        par_kwargs['alpha_min'] = 0
    if 'alpha_max' not in par_kwargs:
        par_kwargs['alpha_max'] = 1

    params = gen_fit_pars(**par_kwargs)

    minimizer = minimize(collective_psd_fit_fun, params,
                         args=(names, freq_dict, model_fun, param_names),
                         kws={'psd_dict': psd_dict,
                              'err_dict': err_dict,
                              'verbose': verbose})

    minimizer.fun = model_fun
    minimizer.param_names = param_names
    minimizer.fun_name = name

    if verbose:
        report_fit(minimizer)

    return minimizer


def is_outlier(red_chi2, N_free, cl=0.95):
    """
    Evaluates the reduced chi² values of the PSD fit and return true if the
    value is considered to be an outlier or false if not.

    Arguments
    ---------
    red_chi2 : array(float) or float
        Numpy array of reduced chi² values. The reduced chi² value is
        calculated as: red_chi² = chi² / N_free
    N_free : array(int) or int
        Number of degrees of freedom.
    conf_level : float
        Confidence level that specifies the range of reduced chi^2 values that
        you are confident about that they come from a normal distribution.

    Returns
    -------
    array of bool

    Note
    ----
    Short
        Say the confidence level is set to 0.95. That means you define every
        datapoint above the 0.975 and below the 0.025 quantile to be an
        outlier. This also means, you expect that, on average, every 20th
        measurement to be an outlier, even if it origins from the expected
        distribution.

    Longer
        In general, the probability that a value lies outside a given
        confidence interval is given by p_out = 1 - conf_level = 2 (1 - p),
        hence p = (1+conf_level)/2. Which, for conf_level = 0.95 means
        p = 0.975.

        The reduced chi2-value (chi2-value divided by the number of degrees of
        freedom) is distributed around a mean = 1; its standard deviation is
        SD = sqrt(2).

        For a large number of degrees of freedom N, the chi2-distribution is a
        normal distribution. The standard error of the mean of normally
        distributed data scales with the square-root of the number of
        measurements, i.e. SEM ~ 1/sqrt(N). Thus, we can ask for the quantile
        for p=0.975 of a normal distribution with mean = 1 and SD = sqrt(2/N)
    """
    # get the quantiles and calculate the upper and lower
    # limit of the reduced chi² values for p = (1+conf_level)/2
    upper_limit = norm.ppf((cl + 1) / 2,
                           loc=1,
                           scale=sqrt(2/N_free)
                           )
    lower_limit = 1 - (upper_limit - 1)
    outlier = (red_chi2 < lower_limit) + (red_chi2 > upper_limit)

    return outlier


class FitResult(object):
    def __init__(self,
                 name=None,
                 D=None,
                 D_err=None,
                 D_unit='V**2/s',
                 f_c=None,
                 f_c_err=None,
                 freq=None,
                 freq_unit='Hz',
                 eval=None,
                 psd_unit='V**2/Hz',
                 residuals=None,
                 chi2=None,
                 redchi2=None,
                 conf_level=None,
                 nfree=None,
                 params=None,
                 minimizer=None,
                 model=''
                 ):
        self.name = name
        self.D = D
        self.D_err = D_err
        self.D_unit = D_unit
        self.f_c = f_c
        self.f_c_err = f_c_err
        self.freq = freq
        self.freq_unit = freq_unit
        self.eval = eval
        self.psd_unit = psd_unit
        self.residuals = residuals
        self.chi2 = chi2
        self.redchi2 = redchi2
        self.conf_level = conf_level
        self.nfree = nfree
        self.params = params
        self.minimizer = minimizer
        self.model = model

    def get_freq_bounds(self):
        return (min(self.freq), max(self.freq))

    def get_D(self, unit='V**2/s'):
        if unit == 'V**2/s':
            conv = 1.0
        else:
            conv = ureg(self.D_unit).to(unit).magnitude
        return self.D * conv

    def get_f_c(self, unit='Hz'):
        if unit == 'Hz':
            conv = 1.0
        else:
            conv = ureg(self.freq_unit).to(unit).magnitude
        return self.f_c * conv

    def is_outlier(self, conf_level=None):
        if not conf_level:
            conf_level = self.conf_level

        return is_outlier(self.redchi2, self.nfree, cl=conf_level)

    def plot_fit(self, axis=None, **kwargs):
        """ plots the fit """
        ax = axis
        if ax is not None:
            fig = ax.figure
        else:
            fig = None

        l = ['plot_all', 'plot_masked', 'plot_errors']
        for key in l:
            try:
                kwargs.pop(key)
            except:
                pass

        if 'color' not in kwargs and self.name in col_dict.keys():
            kwargs['color'] = col_dict[self.name]
        if 'fmt' not in kwargs.keys():
            kwargs['fmt'] = '-'
        if 'linewidth' not in kwargs.keys():
            kwargs['linewidth'] = 1.5
        if 'alpha' not in kwargs.keys():
            kwargs['alpha'] = 0.7
        if 'showLegend' not in kwargs.keys():
            kwargs['showLegend'] = True
        if 'legend_kwargs' not in kwargs:
            lg_kws = {'loc': 3}
            kwargs['legend_kwargs'] = lg_kws
        if 'label' not in kwargs.keys():
            kwargs['label'] = '{} fit to {}'.format(self.model, self.name)
        if 'fontsize' not in kwargs.keys():
            kwargs['fontsize'] = 16

        ax = add_plot_to_figure(fig, self.freq, self.eval,
                                axis=ax, logplot=True, **kwargs)
        ax.grid(which='major')
        plt.setp(ax,
                 xlabel='Frequency (Hz)',
                 ylabel=r'PSD ($\mathsf{V^2/Hz}$)'
                 )
        return ax.figure


class Result(object):
    def __init__(self,
                 dissens,
                 dissens_err,
                 dissens_unit,
                 trap_stiffness,
                 trap_stiffness_err,
                 trap_stiffness_unit,
                 drag,
                 drag_err,
                 drag_unit,
                 excited=None,
                 ):
        self.dissens = dissens
        self.dissens_err = dissens_err
        self.dissens_unit = dissens_unit
        self.trap_stiffness = trap_stiffness
        self.trap_stiffness_err = trap_stiffness_err
        self.trap_stiffness_unit = trap_stiffness_unit
        self.drag = drag
        self.drag_err = drag_err
        self.drag_unit = drag_unit
        self.excited = excited

    def get_dissens(self, unit='nm/mV'):
        if unit != self.dissens_unit:
            conv = ureg(self.dissens_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.dissens * conv

    def get_dissens_err(self, unit='nm/mV'):
        if unit != self.dissens_unit:
            conv = ureg(self.dissens_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.dissens_err * conv

    def get_trap_stiffness(self, unit='pN/nm'):
        if unit != self.trap_stiffness_unit:
            conv = ureg(self.trap_stiffness_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.trap_stiffness * conv

    def get_trap_stiffness_err(self, unit='pN/nm'):
        if unit != self.trap_stiffness_unit:
            conv = ureg(self.trap_stiffness_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.trap_stiffness_err * conv

    def get_drag(self, unit='nN*s/m'):
        if unit != self.drag_unit:
            conv = ureg(self.drag_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.drag * conv

    def get_drag_err(self, unit='nN*s/m'):
        if unit != self.drag_unit:
            conv = ureg(self.drag_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.drag_err * conv

    def get_force_factor(self, unit='pN/mV'):
        """TODO"""
        beta = self.get_dissens(unit='nm/mV')
        kappa = self.get_trap_stiffness(unit='pN/nm')
        if unit != 'pN/mV':
            conv = ureg('pN/mV').to(unit).magnitude
        else:
            conv = 1.0
        alpha = beta * kappa * conv
        return alpha

    def get_force_factor_err(self, unit='pN/mV'):
        """ TODO """
        beta = self.get_dissens(unit='nm/mV')
        kappa = self.get_trap_stiffness(unit='pN/nm')
        beta_err = self.get_dissens_err(unit='nm/mV')
        kappa_err = self.get_trap_stiffness_err(unit='pN/nm')
        alpha = self.get_force_factor(unit=unit)

        alpha_err = (beta_err/beta + kappa_err/kappa) * alpha
        return alpha_err


class PSDFit(object):
    def __init__(self,
                 psdm,
                 bounds=None,
                 conf_level=0.95,
                 verbose=False,
                 **kwargs):
        """
        PSDFit object manages the several fitting procedures of power spectral
        densities of one measurement.

        When initialized, the *analytical_lorentzian_fit* method is called,
        which tries to find optimal bounds and initial values for the fits.

        Arguments
        ---------
        psdm : PSDMeasurement
            Object holding the power spectral densities of x, y and z.
        bounds : (fmin, fmax) or dict {name: (fmin, fmax), ...}
            defines the lower and upper frequency bound. All data outside the
            interval are masked. A dict specifies the bound for particular
            axes. The tuple will put the same boundary to each spectrum.
        kwargs
            Keyword arguments handed over to **analytical_lorentzian_fit()**.
        """
        self.psdm = psdm

        self.active_calibration = psdm.active_calibration
        self.names = psdm.names

        self.anal_fits = {}

        self.fit_kwargs = {}
        self.fits = {}
        
        self.conf_level = conf_level

        if bounds is not None:
            self.set_bounds(bounds)
        self.analytical_lorentzian_fit(verbose=verbose, **kwargs)

        self.pc_results = {}
        self.ac_results = {}

        self.exclude_freq_outside = self.set_bounds

    def __getattr__(self, name):
        """
        Shortcut for PSDM.attribute
        """
        if hasattr(self.psdm, name):
            return getattr(self.psdm, name)
        else:
            raise AttributeError(name)

    @property
    def freq_x(self):
        return self.get_freq('x')

    @property
    def freq_y(self):
        return self.get_freq('y')

    @property
    def freq_z(self):
        return self.get_freq('z')

    def get_freq(self, name, **kwargs):
        """
        Return the frequency vector of the psd with 'name'.

        Arguments
        ---------
        name : str
            name of the axis.
        kwargs can be:
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        """
        return self.psdm.get_freq(name=name, **kwargs)

    @property
    def psd_x(self):
        return self.get_psd('x')

    @property
    def psd_y(self):
        return self.get_psd('y')

    @property
    def psd_z(self):
        return self.get_psd('z')

    def get_psd(self, name, **kwargs):
        """
        Return the psd vector of the psd with 'name'.

        Arguments
        ---------
        name : str
            name of the axis.
        kwargs can be:
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        """
        return self.psdm.get_psd(name=name, **kwargs)

    @property
    def err_x(self):
        return self.get_psd_err('x')

    @property
    def err_y(self):
        return self.get_psd_err('y')

    @property
    def err_z(self):
        return self.get_psd_err('z')

    def get_psd_err(self, name, **kwargs):
        """
        Return the error vector of the psd with 'name'.

        Arguments
        ---------
        name : str
            name of the axis.
        kwargs can be:
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        """
        return self.psdm.get_psd_err(name=name, **kwargs)

    def set_bounds(self, bounds, names=None, reset_mask=False):
        """
        Set the frequency boundaries (fmin, fmax) at which to consider the psd
        values.

        Arguments
        ---------
        bounds : tuple (fmin, fmax) or dict {'name': (fmin, fmax), ...}
            Defines the lower and upper frequency bound. All data outside the
            interval is masked.
        names : list of strings
            Names to apply the tuple bounds to.
        reset : bool
            Resets the mask of the psd.
        """
        if hasattr(bounds, 'keys'):
            for name in bounds:
                fmin, fmax = bounds[name]
                self.psdm.exclude_freq_outside(fmin, fmax, names=name,
                                               reset_mask=reset_mask)
        else:  # either use bounds for all names or for 'names'
            if names and not isinstance(names, list):
                names = [names]

            for name in names if names else self.names:
                fmin, fmax = bounds
                self.psdm.exclude_freq_outside(fmin, fmax, names=name,
                                               reset_mask=reset_mask)

    def exclude_freq(self, f_exclude, names=None):
        """
        Exclude frequencies.

        Arguments
        ---------
        f_exclude : float or list of floats or dictionary that specifies the
            frequencies for each name. 'names' then has no effect.
        names : list of str
            Names to apply the exclusion of the list or float 'f_exclude'.
        """
        if hasattr(f_exclude, 'keys'):
            for name in f_exclude:
                self.psdm.exclude_freq(f_exclude[name], names=name)
        else:
            if names and not isinstance(names, list):
                names = [names]

            self.psdm.exclude_freq(f_exclude, names=names)

    def analytical_lorentzian_fit(self, names=None, verbose=False, **kwargs):
        """
        Calculates the analytical solution to the least squares problem of a
        simple Lorentzian.

        This function calls *calc_anal_lsq_opt()* to find a frequency interval
        where the data fits a lorentzian well.

        Arguments
        ---------
        names : [str, ...] or 'all'
            Defines for which PSD(s) the analytical lorentzian is calculated.
            Default is 'all'.
        verbose : bool
            Give more info.

        Keyword Arguments
        -----------------
        Keyword arguments passed over to *calc_anal_lsq_opt()*.

        References
        ----------
        K. Berg-Sørensen and H. Flyvbjerg
            Power spectrum analysis for optical tweezers
            Rev. Sci. Instr., vol. 75, 2004
        S.F. Nørrelykke and H. Flyvbjerg
            Power spectrum analysis with least-squares fitting: Amplitude bias
            and its elimination, with application to optical tweezers and
            atomic force microscope cantilevers
            Rev. Sci. Instr.,vol. 81, 2010

        See Also
        --------
        calc_anal_lsq_opt
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            # calculate de-biasing factor
            n = self.psdm.psds[name].N_avg
            debias = n / (n + 1)

            if self.active_calibration and name is self.psdm.ex_axis:
                f_ex = self.psdm.get_ex_freq(unit='Hz')
                self.psdm.exclude_freq(f_ex, names=self.psdm.ex_axis)

            freq = self.get_freq(name, unit='Hz')
            psd = self.get_psd(name, unit='V**2/Hz')
            (f_c, D, chi2, bounds) = calc_anal_lsq_opt(freq, psd,
                                                       verbose=verbose,
                                                       **kwargs)
            # check for nans at the output of calc_anal_lsq_opt
            N_nans = sum(1 * isnan(val) for val in [f_c, D, chi2])
            if N_nans > 0:
                model = None
                f_c = 1000
                D = 1
                chi2 = 0
                bounds = (min(freq), max(freq))
                if verbose:
                    print('Analytical fit did not succeed. Using the '
                          'following default values: '
                          'f_c=1000 hz, D=1.0 V²/s, chi2=0, bounds={}'
                          ''.format(bounds))
            else:
                model = 'analytical least squares'

            keep = logical_and(bounds[0] <= freq, freq <= bounds[1])
            freq_keep = freq[keep]
            N_free = len(freq_keep) - 2
            lf_eval = lorentzian_psd(freq_keep, D * debias, f_c)

            FR = FitResult(name=name,
                           D=D * debias,
                           D_err=None,
                           D_unit='V**2/s',
                           f_c=f_c,
                           f_c_err=None,
                           freq=freq_keep,
                           freq_unit='Hz',
                           eval=lf_eval,
                           psd_unit='V**2/Hz',
                           residuals=None,
                           chi2=chi2,
                           redchi2=chi2/N_free,
                           nfree=N_free,
                           params=None,
                           minimizer=None,
                           model=model)

            self.anal_fits.update({name: FR})

    def plot_anal_fits(self, names=None, axis=None, **kwargs):
        """
        Plot the analytical least squares fits.

        Arguments
        ---------
        names : str
            specify which axis to plot. E.g. names='xy', would plot only x and
            y.
        axis : matplotlib:Axis
            Specify in which axis to add the plot
        **kwargs
            Keyword arguments handed over to *plot_fit()*.
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            fig = self.anal_fits[name].plot_fit(axis=axis, **kwargs)
            if axis is None:
                axis = fig.axes[0]

        return fig

    def collective_psd_fit(self,
                           names=None,
                           model='lorentzian',
                           fitreport=False,
                           aliasing=False,
                           N_alias=9,
                           verbose=False,
                           **kwargs):
        """
        Fits the power spectra to the selected model function collectively,
        which is modified by a low-pass filter. The parameters alpha and
        f3dB of the low-pass filter are the same for all psds.
        This can be used to find alpha and f3dB.

        Arguments
        ---------
        names : str or list of str
            names of the PSDs to include in the collective fit, defaults to all
            available names.
        model : {'lorentzian', 'hydro'}
            Define the model function to be used.
        aliasing : bool
            Whether to include aliasing. If True, the sample frequency of the
            psdmeasurment is used.
        N_alias : int
            Number that defines how many ranges of f_sample should be taken
            into account. The default N=9 give a very good approximation with
            deviations of less than 0.5%
        verbose : bool
            Give more info.

        Keyword Arguments
        -----------------
        Keyword arguments are passed to make_collective_psd_fit. They can be
        initial guess values for the parameters f3dB and alpha and D and f_c
        where the axes for the latter two are specified by '_name',
        i.e. D_x=1 would specify the initial guess for the diffusion constant
        of the axis 'x'. Further, the suffixes '_min' and '_max' identify lower
        and upper boundaries for the fitting.
        """
        if model.lower() is 'hydro':
            expset = self.psdm.exp_setting
            h_kw = {'radius': expset.get_radius(unit='m'),
                    'height': expset.get_height(unit='m'),
                    'temp': expset.get_temp(unit='K'),
                    'rho': expset.get_density_particle(unit='kg/m**3'),
                    'density_med': expset.get_density_medium(unit='kg/m**3'),
                    'viscosity': expset.get_viscosity(unit='Pa*s')}
            kwargs.update(h_kw)

        freq_dict = OrderedDict()
        psd_dict = OrderedDict()
        err_dict = OrderedDict()

        fs_ = []

        if names and not isinstance(names, list):
            names = [names]

        names = names if names else self.names

        for name in names:
            freq_dict[name] = self.get_freq(name, unit='Hz')
            psd_dict[name] = self.get_psd(name, unit='V**2/Hz')
            err_dict[name] = self.get_psd_err(name, unit='V**2/Hz')
            if 'D_' + name not in kwargs:
                kwargs['D_' + name] = self.anal_fits[name].get_D(unit='V**2/s')
            if 'f_c_' + name not in kwargs:
                kwargs['f_c_' + name] = self.anal_fits[name].get_f_c(unit='Hz')
            fs_.append(self.psdm.psds[name].get_f_sample(unit='Hz'))

        if len(set(fs_)) > 1:
            warnings.warn('Modelling PSDs with different sampling rates '
                          'is not yet supported! Modelling aliasing will not '
                          'work correctly.')
        minimizer = make_collective_psd_fit(freq_dict,
                                            psd_dict,
                                            err_dict,
                                            model=model,
                                            aliasing=aliasing,
                                            f_sample=fs_[0],
                                            N_alias=9,
                                            verbose=verbose,
                                            **kwargs)
        params = minimizer.params

        model_ = collective_psd_fit_fun(params,
                                        names,
                                        freq_dict,
                                        minimizer.fun,
                                        minimizer.param_names)
        start = 0
        npara_per_axis = 2  # f_c and D
        ntotal = len(minimizer.residual)

        for name in names:
            freq = self.get_freq(name)
            len_ = len(freq)
            mod_ = model_[start:start+len_]
            res_ = minimizer.residual[start:start+len_]
            start += len_

            chi2 = (res_**2).sum()
            ndata = len_
            # f3dB and alpha are shared among all datapoints
            nshared = 2 * ndata / ntotal
            nfree = ndata - npara_per_axis - nshared
            redchi2 = chi2 / nfree

            fres = FitResult(name=name + ' (cltv)',
                             D=params['D_' + name].value,
                             D_err=params['D_' + name].stderr,
                             f_c=params['f_c_' + name].value,
                             f_c_err=params['f_c_' + name].stderr,
                             freq=freq,
                             eval=mod_,
                             residuals=res_,
                             chi2=chi2,
                             redchi2=redchi2,
                             conf_level=self.conf_level,
                             nfree=nfree,
                             params=params,
                             minimizer=minimizer,
                             model='collective lp-filtered ' + model +' psd fit'
                             )

            self.fits.update({name: fres})
            
        if fitreport and not verbose:
            report_fit(minimizer, show_correl=False)

    def setup_fit(self,
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
                  **fit_psd_kws):
        """
        Setup the psd fitting parameters.

        Arguments
        ---------
        names : str or list of str
            Defines the name(s) of the fits to be set up.
        model : {'lorentzian', 'hydro'}
            Define the model function to be used.
        lp_filter : bool
            If True the selected model will be modified by a low-pass
            filter.
        lp_fixed : bool
            If true the low pass filter parameters alpha and f3dB
            will be fixed. Expects the keywords 'alpha' and 'f3dB' and their
            respective values as kwargs. If they are not provided they are set
            to alpha = 0.5 and f3dB=10e3 Hz.
        aliasing : bool
            If aliasing should be considered or not.
        f_sample : float
            If aliasing is True, the sampling frequency of the original signal
            is needed to model the aliasing. If None, it is set to 2 * f_max,
            where f_max is the highest frequency in the frequency vector.
        N_alias : int
            Number that defines how many ranges of f_sample should be taken
            into account. The default N=9 give a very good approximation with
            deviations of less than 0.5%
        debias : bool
            Whether to account for biasing after the fitting. This only effects
            the Diffusion coefficient by a factor of N_avg / (N_avg + 1).
            N_avg is taken from the PSDMeasurement object.
        f3dB : float
            Cut-off frequency of the low-pass filter. Either the initial guess
            or the fixed value if lp_fixed is True.
        alpha : float
            Filter efficiency. Only a fraction (0 <= alpha <= 1) of the signal
            that is not low pass filtered. Thus, alpha = 1, will produce no
            low-pass filter at all, whereas alpha = 0, will produce a 1st-order
            low-pass filter.
        **fit_psd_kws : keyword arguments passed to
            pyotc.psd_fitting.fit_psd().
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            kws = fit_psd_kws
            kws['model'] = model.lower()
            kws['lp_filter'] = lp_filter
            kws['lp_fixed'] = lp_fixed
            kws['aliasing'] = aliasing
            if aliasing and f_sample is None:
                try:
                    f_sample = self.psdm.psds[name].get_f_sample(unit='Hz')
                except:
                    pass
            kws['f_sample'] = f_sample
            kws['N_alias'] = N_alias
            kws['debias'] = debias
            if f3dB:
                kws['f3dB'] = f3dB
            if alpha:
                kws['alpha'] = alpha
            self.fit_kwargs[name] = kws

    def fit_psd(self,
                name,
                dynamic_bounds=False,
                bounds=None,
                f_exclude=None,
                use_heights=False,
                plot_fit=False,
                plot_axis=None,
                plot_masked_data=True,
                fitreport=False,
                verbose=False,
                **kwargs):
        """
        Fit psd of data 'name' to the selected model.

        Consider setting up the fit first by self.setup_fit().

        Arguments
        ---------
        name : str
            name of the axis to be fitted.
        dynamic_bounds : bool
            Whether to determine the frequency bounds from the analytical fit.
            If True the upper bound is set one order of magnitude
            higher than the analytically determined corner freuqency. The
            lower bound matches the optimal minimum frequency from the
            calc_anal_lsq_opt() function, if this frequency is less than one
            order of magnitude lower than the analytically determined corner
            freuqency.
        bounds : tuple (fmin, fmax) as limits for fitting
        f_exclude : float or list
            frequency or a list of frequencies that should be excluded from
            the fit.
        use_heights : bool
            Only important for hydrodyn. cor. fitting of a psd. If True, the
            internally set height is used otherwise it is set to infinity.
        plot_fit : bool
            Whether to plot the fit.
        plot_axis : matplotlib.Axis or None
        plot_masked : bool
            Whether to plot masked data as well.
        fitreport : bool
            Whether to print a report of the fit or not.
        verbose : bool
            Print details for debugging.

        Keyword Arguments
        -----------------
        height : float
            If model is 'hydro'
            Distance between the bead center and the surface,
            if it is not specified infinity is used.
        density_med : float
            Density of the medium in kg/m³. Only used if model is 'hydro'. If
            None, the density of water is calculated at the temperature
            specified in the measurement object.
        viscosity : float
            Viscosity of the medium. Only used if model is 'hydro'. If
            None, the viscosity of water is calculated at the temperature
            specified in the measurement object.

        Note
        ----
        Other Keyword arguments can also be initial guesses of 'D' and 'f_c',
        one can also set their limits by specifying '_min' and '_max' for the
        respective fitting parameter (see also: gen_fit_pars()).

        If model is 'hydro', the values of the particle radius, temperature and
        rho (the density of the particle) are taken from the PSDMeasurement-
        object.

        Setting bounds for a certain 'name' resets the former mask of the
        corresponding PSD-object. All other masked values at specific
        frequencies must be defined again via f_exclude.
        """
        if bounds or f_exclude:
            self.psdm.psds[name].reset_mask()

        if dynamic_bounds:
            anal_fit = self.anal_fits[name]
            f_c = anal_fit.f_c
            opt_fmin = anal_fit.get_freq_bounds()[0]
            # take f_min from analytical fit only if this bound is one
            # order of magnitude smaller than the corner frequency
            fmin = opt_fmin if opt_fmin < f_c / 10 else f_c / 10
            bounds = (fmin, 10 * f_c)

        if bounds:
            fmin, fmax = bounds
            self.psdm.exclude_freq_outside(fmin, fmax, names=name)
        if f_exclude:
            self.psdm.exclude_freq(f_exclude, names=name)

        if self.active_calibration and name == self.psdm.ex_axis:
            f_ex = self.psdm.get_ex_freq(unit='Hz')
            self.psdm.exclude_freq(f_ex, names=name)
            isexcited = True
        else:
            isexcited = False

        # create fit options if not done yet
        if name not in self.fit_kwargs:
            self.setup_fit(names=name)
            warnings.warn('Fit was not setup. Using default settings. '
                          'See setup_fit()')
        fitkwargs = self.fit_kwargs[name]

        kwargs.update(fitkwargs)

        # should the fit be de-biased?
        N_avg = self.psdm.psds[name].N_avg
        if kwargs.pop('debias'):
            deb = N_avg / (N_avg + 1)
        else:
            deb = 1

        # set initial values for D and f_c if not provided
        if 'D' not in kwargs:
            kwargs['D'] = self.anal_fits[name].D
        if 'f_c' not in kwargs:
            kwargs['f_c'] = self.anal_fits[name].f_c

        if kwargs['model'] == 'hydro':
            expset = self.psdm.exp_setting

            radius = expset.get_radius(unit='m')
            temp = expset.get_temp(unit='K')
            density_p = expset.get_density_particle(unit='kg/m**3')
            density_m = expset.get_density_medium(unit='kg/m**3')
            viscosity = expset.get_viscosity(unit=('Pa*s'))
            lateral = self.psdm.psds[name].is_lateral()

            if use_heights:
                height = kwargs['height'] = expset.get_height(unit='m')
            else:
                height = inf

            kwargs['radius'] = radius
            kwargs['temp'] = temp
            kwargs['rho'] = density_p
            kwargs['density_med'] = density_m
            kwargs['viscosity'] = viscosity
            kwargs['lateral'] = lateral

            s_pars = ('Radius = {0:1.3e} m\t height = {1:1.3e} m \t '
                      'T = {2:1.3f} K\t density = {3:1.3f}\n'
                      'density medium = {4:1.3f} kg/m**3 \t '
                      'viscosity = {5:1.3f} Pa*s'
                      ''.format(radius, height, temp, density_p,
                                density_m, viscosity))
        else:
            s_pars = ''

        lp_filter = kwargs['lp_filter']
        lp_fixed = kwargs['lp_fixed']

        if fitreport or verbose:
            print('----------')
            print(" '{0:s}' ".format(name).center(10))
            print('----------')

        freq = self.get_freq(name, unit='Hz')
        psd = self.get_psd(name, unit='V**2/Hz')
        err = self.get_psd_err(name, unit='V**2/Hz')

        modelfit = fit_psd(freq,
                           psd,
                           psd_err=err,
                           N_avg=N_avg,
                           verbose=verbose,
                           **kwargs)

        modelname = modelfit.model.name.replace('Model(', '').replace(')', '')
        if fitreport or verbose:
            print('Model: {0:s}'.format(modelname, name))
            print('Flags and parameters:')
            print('Low pass filter: {}\t fixed: {}'
                  ''.format(lp_filter, lp_fixed))
            if lp_filter:
                print('f_3dB = {0:1.1f}\t alpha = {1:1.3f}'
                      ''.format(kwargs['f3dB'], kwargs['alpha']))
            if s_pars:
                print(s_pars)
            print('De-bias factor = {0:1.3f}'.format(deb))

            if not verbose:
                report_fit(modelfit)

        params = modelfit.params
        FR = FitResult(name=name,
                       D=params['D'].value * deb,
                       D_err=params['D'].stderr * deb,
                       D_unit='V**2/s',
                       f_c=params['f_c'].value,
                       f_c_err=params['f_c'].stderr,
                       freq=freq,
                       freq_unit='Hz',
                       eval=1 / modelfit.eval() * deb,
                       psd_unit='V**2/Hz',
                       residuals=modelfit.residual,
                       chi2=modelfit.chisqr,
                       redchi2=modelfit.redchi,
                       conf_level=self.conf_level,
                       nfree=modelfit.nfree,
                       params=params,
                       minimizer=modelfit,
                       model=modelname
                       )
        self.fits[name] = FR

        # The result of excited axis can be accessed without knowing the name
        if isexcited:
            if 'excited_axis' in self.fits.keys():
                warnings.warn('Either the fit of the excited axis was repeated'
                              ' or there are more than one excited axis. '
                              'The "excited_axis" key is overwritten.'
                              'Calling the fit intrisically by "excited_axis" '
                              'and not by the name of the axis can lead to '
                              'unwanted results.')
            self.fits['excited_axis'] = FR

        if plot_fit:
            self.plot_fits(names=name,
                           plot_axis=plot_axis,
                           plot_masked_data=plot_masked_data)

    def fit_psds(self,
                 names=None,
                 dynamic_bounds=False,
                 bounds=None,
                 f_exclude=None,
                 use_heights=False,
                 calc_results=True,
                 plot_fits=False,
                 plot_axis=None,
                 fitreport=False,
                 verbose=False,
                 **kwargs):
        """
        Fit psds.

        Same method as fit_psd. this one iterates through the axes 'names' or
        all available axes in self.names.

        Arguments
        ---------
        names : str or list of str
            names to run the fiiting routine all available once are used if
            None.
        dynamic_bounds : bool
            Whether to determine the frequency bounds from the analytical fit.
            If True the upper bound is set one order of magnitude
            higher than the analytically determined corner freuqency. The
            lower bound matches the optimal minimum frequency from the
            calc_anal_lsq_opt() function, if this frequency is less than one
            order of magnitude lower than the analytically determined corner
            freuqency.
        bounds : dict
            Dictionary with names of the axes as keys and tuples of fmin and
            fmax as values.
        f_exclude : dict
            Dictionary with names of the axes as keys and lists of frequencies
            to be excluded.
        use_heights : bool
            Only important for hydrodyn. cor. fitting of a psd. If True, the
            internally set height is used otherwise it is set to infinity.
        calc_results : bool
            Calculate the pc and ac results if True.
        plot_fits : bool
            Whether to plot the fits.
        plot_axis : matplotlib.Axis or None
        fitreport : bool
            Whether to print a report of the fit or not.
        verbose : bool
            Print details for debugging.
        """
        if fitreport or verbose:
            print('========================================================')

        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            if hasattr(bounds, 'keys'):
                try:
                    bnds = bounds[name]
                except:
                    bnds = None
            else:
                bnds = bounds

            if hasattr(f_exclude, 'keys'):
                try:
                    f_ex = f_exclude[name]
                except:
                    f_ex = None
            else:
                f_ex = f_exclude

            self.fit_psd(name=name,
                         dynamic_bounds=dynamic_bounds,
                         bounds=bnds,
                         f_exclude=f_ex,
                         use_heights=use_heights,
                         fitreport=fitreport,
                         verbose=verbose,
                         **kwargs)
            if calc_results:
                self.calc_pc_results(names=name, verbose=verbose)

        if plot_fits:
            self.plot_fits(names=names,
                           plot_axis=plot_axis)

        if fitreport or verbose:
            print('========================================================')

        if self.active_calibration:
            if all((name in self.fits) for name in self.names):
                self.calc_ac_results(verbose=verbose)

        if verbose:
            self.print_pc_results()
            self.print_ac_results()

        if plot_fits:
            plt.show()

    def plot_fits(self,
                  names=None,
                  plot_axis=None,
                  plot_data=True,
                  plot_masked_data=True,
                  save_plots=False,
                  save_as='png',
                  directory=None,
                  filename=None,
                  **kwargs):
        """
        Plots the fits.
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            if plot_data:
                fig = self.psdm.plot_psds(names=name,
                                          axis=plot_axis,
                                          **kwargs)

                if plot_axis is None:
                    plot_axis = fig.axes[0]

                if plot_masked_data:
                    self.psdm.plot_psds(names=name,
                                        axis=plot_axis,
                                        plot_masked=True,
                                        **kwargs)

            fig = self.fits[name].plot_fit(axis=plot_axis, **kwargs)
            if plot_axis is None:
                plot_axis = fig.axes[0]

        if save_plots:
            directory = directory if directory else './'
            ptf = ''.join([directory, filename, '.', save_as])
            fig.patch.set_alpha(1.0)
            fig.patch.set_color('white')
            fig.savefig(ptf, dpi=150, format=save_as)

        return fig

    def is_outlier(self, names=None, conf_level=None):
        """
        Do the fits lie within the given confidence interval?
        """
        if not conf_level:
            conf_level = self.conf_level

        if names and not isinstance(names, list):
            names = [names]

        names = names if names else self.names

        return {name: self.fits[name].is_outlier(conf_level)
                for name in names}

    def get_basepower(self, unit='V**2'):
        """
        Determines the basepower and its error from the fit.
        The error is calculated and should be understood as an upper error
        estimate. A correct calculation should consider the covariances
        from the fitting parameters.
        As far as the excitation amplitude is at least two orders of magnitude
        higher as the basepower, the impact that this error has is of minor
        interest.
        """
        if unit != self._ex_power_unit:
            conv = ureg(self._ex_power_unit).to(unit).magnitude
        else:
            conv = 1.0

        ex_axis = self.psdm.ex_axis
        fit = self.fits[ex_axis]
        f_ex = self.psdm.get_ex_freq(unit='Hz')
        df = self.psdm.psds[ex_axis].df

        # evaluate the fit at the excitation frequency
        # take inverse because fit is done inversely
        bp = (1 / fit.minimizer.eval(f=array([f_ex]))[0]) * df * conv
        bp_err = bp / sqrt(self.psdm.psds[ex_axis].N_avg) * conv

        return (bp, bp_err)

    def calc_pc_results(self, names=None, verbose=False):
        """
        Calculate the passive calibration results.

        The following values are calculated:

         - displacement sensitivity
         - trap stiffness
         - Stokes drag

        the calculations assume Stokes' drag of a sphere of the given radius,
        viscosity and temperature.
        """
        T = self.psdm.exp_setting.get_temp(unit='K')
        dT = self.psdm.exp_setting.get_temp_err(unit='K')

        drag = self.psdm.get_stokes_drag(unit='N*s/m')
        ddrag = self.psdm.get_stokes_drag_err(unit='N*s/m')

        D_0 = (k_B * T / drag)  # Nm/K * K / (Ns/m) = m²/s
        dD_0 = (dT / T + ddrag / drag) * D_0

        if verbose:
            print('passive calibration results:')
            print('T = {0:1.3e} +/- {1:1.3e} K'.format(T, dT))
            print('Stokes drag = {0:1.3e} +/- {1:1.3e} Ns/m'
                  ''.format(drag, ddrag))
            print('D_0 = k_B * T / (6*pi*eta*R) = {0:1.3e} +/- {1:1.3e} m^2/s'
                  ''.format(D_0, dD_0))

        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            if name not in self.fits.keys():
                raise Exception("Fit for '{}' not present."
                                "Try to fit the PSD first.".format(name))
            D_V = self.fits[name].D  # V**2/s
            dD_V = self.fits[name].D_err

            f_c = self.fits[name].f_c  # Hz
            df_c = self.fits[name].f_c_err

            # dissens
            beta_pc = (D_0 / D_V)**0.5 * 1e6  # m/V = 1e6 nm/mV
            dbeta_pc = (0.5 * dD_0 / D_0 + dD_V / D_V) * beta_pc

            # trap_stiffness
            # [kappa] = 1/s * Ns/m = N/m = 1e3 pN/nm
            kappa_pc = (2 * pi * f_c * drag) * 1e3
            dkappa_pc = (df_c / f_c + ddrag / drag) * kappa_pc

            pc_result = Result(dissens=beta_pc,
                               dissens_err=dbeta_pc,
                               dissens_unit='nm/mV',
                               trap_stiffness=kappa_pc,
                               trap_stiffness_err=dkappa_pc,
                               trap_stiffness_unit='pN/nm',
                               drag=drag * 1e9,  # Ns/m = 1e9 nNs/m
                               drag_err=ddrag * 1e9,
                               drag_unit='nN*s/m',
                               excited=False,
                               )
            self.pc_results.update({name: pc_result})

            if verbose:
                print(name)
                print('D_V = {0:1.3e} +/- {1:1.3e} V^2/s'.format(D_V, dD_V))
                print('f_c = {0:1.3e} +/- {1:1.3e} Hz'.format(f_c, df_c))
                print('beta_pc = {0:1.3e} +/- {1:1.3e} nm/mV'
                      ''.format(beta_pc, dbeta_pc))
                print('kappa_pc = {0:1.3e} +/- {1:1.3e} pN/nm'
                      ''.format(kappa_pc, dkappa_pc))

        if self.active_calibration and self.psdm.ex_axis in self.pc_results:
            self.pc_results['excited_axis'] = self.pc_results[self.psdm.ex_axis]

    def calc_ac_results(self, verbose=False):
        """
        Calculates the values:
         - displacement sensitivity
         - trap stiffness
         - drag
        for the given PSDM.
        """
        names = self.names.copy()
        for name in names:
            if name not in self.fits.keys():
                raise Exception("Fit for '{}' not present."
                                "Try to fit the PSD first.".format(name))
        # name of ac-excited PSD
        ex_axis = self.psdm.ex_axis
        # TODO axial ac-calibration
        if self.psdm.psds[ex_axis].direction is 'axial':
            warnings.warn('Axial excitation involves height dependent-changes '
                          'of the drag coefficient that are dependent on the '
                          'amplitude of excitation itself. These '
                          'dependencies are not yet accounted for. Thus '
                          'it is recommended to excite in the lateral '
                          'direction. Ignoring this you accept the risk of '
                          'large errors. In other words, be sure you really '
                          'know what you do.')
        # excited frequency
        f_ex = self.psdm.ex_freq
        # excitation amplitude in meters
        A = self.psdm.get_ex_amplitude(unit='m')
        dA = self.psdm.get_ex_amplitude_err(unit='m')
        # power at excited frequency in V**2/Hz
        P = self.psdm.get_ex_power(unit='V**2')
        dP = self.psdm.get_ex_power_err(unit='V**2')
        # thermal (base) power density at excited frequency
        base, dbase = self.get_basepower(unit='V**2')
        # corner frequency from fit to excited axis
        f_c = self.fits[ex_axis].f_c
        df_c = self.fits[ex_axis].f_c_err
        # theoretical power for that excited frequency in m²
        W_th = A**2 / (2 * (1 + (f_c/f_ex)**2))
        dW_th = (2 * (dA / A + df_c / f_c)) * W_th
        # actual (measured) power in V²
        W_exp = (P - base)
        dW_exp = (dP + dbase)

        if verbose:
            print('ex amp: {0:1.3e} +/- {1:1.3e} m'.format(A, dA))
            print('ex power: {0:1.3e} +/- {1:1.3e} V**2'
                  ''.format(P, dP))
            print('base power: '
                  '{0:1.3e} +/- {1:1.3e} V**2'.format(base, dbase))
            print('W_th: {0:1.3e} +/- {1:1.3e} m**2'.format(W_th, dW_th))
            print('W_exp: {0:1.3e} +/- {1:1.3e} V**2'.format(W_exp, dW_exp))

        # displacement sensitivity for the excited axis
        beta_ac = (W_th / W_exp)**0.5  # sqrt(m²/V²) = m/V
        dbeta_ac = (0.5 * (dW_th / W_th + dW_exp / W_exp)) * beta_ac
        # Temperature in Kelvin
        T = self.psdm.exp_setting.get_temp(unit='K')
        dT = self.psdm.exp_setting.get_temp_err(unit='K')
        # diffusion constant determined from fit (V²/s)
        D_ac = self.fits[ex_axis].D
        dD_ac = self.fits[ex_axis].D_err
        # drag in physical units Nm/K * K / (m²/V² * V²/s) = Ns/m
        gamma_ac = k_B * T / (beta_ac**2 * D_ac)
        dgamma_ac = (dT / T + 2 * dbeta_ac / beta_ac + dD_ac / D_ac) * gamma_ac
        # trap stiffness in physical units (N/m)
        kappa_ac = 2 * pi * f_c * gamma_ac
        dkappa_ac = (df_c / f_c + dgamma_ac / gamma_ac) * kappa_ac

        if verbose:
            print('beta_ac = {0:1.3e} +/- {1:1.3e} nm/mV'
                  ''.format(beta_ac * 1e6, dbeta_ac * 1e6))
            print('T = {0:1.3e} +/- {1:1.3e} K'.format(T, dT))
            print('D_ac = {0:1.3e} +/- {1:1.3e} V**2/s'.format(D_ac, dD_ac))
            print('gamma_ac = {0:1.3e} +/- {1:1.3e} nN*s/m'
                  ''.format(gamma_ac * 1e9, dgamma_ac * 1e9))
            print('kappa_ac = {0:1.3e} +/- {1:1.3e} pN/nm'
                  ''.format(kappa_ac * 1e3, dkappa_ac * 1e3))

        ac_result = Result(dissens=beta_ac * 1e6,
                           dissens_err=dbeta_ac * 1e6,
                           dissens_unit='nm/mV',
                           drag=gamma_ac * 1e9,
                           drag_err=dgamma_ac * 1e9,
                           drag_unit='nN*s/m',
                           trap_stiffness=kappa_ac * 1e3,
                           trap_stiffness_err=dkappa_ac * 1e3,
                           trap_stiffness_unit='pN/nm',
                           excited=True
                           )
        self.ac_results[ex_axis] = ac_result
        self.ac_results['excited_axis'] = ac_result

        names.remove(ex_axis)

        for name in names:
            D = self.fits[name].D  # V**2/s
            dD = self.fits[name].D_err
            f_c = self.fits[name].f_c  # Hz
            df_c = self.fits[name].f_c_err

            beta = beta_ac * sqrt(D_ac / D)  # m/V
            dbeta = (dbeta_ac / beta_ac + 0.5 * (dD_ac / D_ac + dD / D)) * beta

            kappa = 2 * pi * f_c * gamma_ac  # N/m
            dkappa = (df_c / f_c + dgamma_ac / gamma_ac) * kappa

            ac_result = Result(dissens=beta * 1e6,
                               dissens_err=dbeta * 1e6,
                               dissens_unit='nm/mV',
                               drag=gamma_ac * 1e9,
                               drag_err=dgamma_ac * 1e9,
                               drag_unit='nN*s/m',
                               trap_stiffness=kappa * 1e3,
                               trap_stiffness_err=dkappa * 1e3,
                               trap_stiffness_unit='pN/nm',
                               excited=False
                               )
            self.ac_results[name] = ac_result

    def print_ac_results(self):
        self.print_results(active_calibration=True)

    def print_pc_results(self):
        self.print_results(active_calibration=False)

    def print_results(self, active_calibration=True):
        """
        Print the results of the psd fits.
        """
        if active_calibration:
            if not self.active_calibration:
                print("PSD fit is set to 'passive calibration'. Check the "
                      "measurement-object's active_calibration attribute and "
                      "the ones starting with ex_.")
                return None
            s = 'ac_'
            s2 = 'Active calibration results'
        else:
            s = 'pc_'
            s2 = 'Passive calibration results'

        print('========================================================')
        print(s2)

        for name in self.names:
            # get the results of the passive or active calibration calculations
            result = getattr(self, s + 'results')[name]
            print('--------------------------------------------------------')
            print("'{}'".format(name))
            print('--------------------------------------------------------')
            print('Displacement sensitivity: {0:1.3e} +/- {1:1.3e} {2}'
                  ''.format(result.get_dissens(unit='nm/mV'),
                            result.get_dissens_err(unit='nm/mV'),
                            result.dissens_unit))
            print('Trap stiffness: {0:1.3e} +/- {1:1.3e} {2}'
                  ''.format(result.get_trap_stiffness(unit='pN/nm'),
                            result.get_trap_stiffness_err(unit='pN/nm'),
                            result.trap_stiffness_unit)
                  )
            print('Drag: {0:1.3e} +/- {1:1.3e} {2}'
                  ''.format(result.get_drag(unit='nN*s/m'),
                            result.get_drag_err(unit='nN*s/m'),
                            result.drag_unit))
        print('========================================================')

    def write_results_to_file(self, directory=None, fname=None, fext='.txt'):
        """
        Write fit results to a file.
        """
        if directory is None:
            try:
                directory = self.psdm.directory
            except:
                directory = './'
        if fname is None:
            try:
                fname = self.psdm.paramfile
            except:
                dstr = time.strftime("%Y-%m-%d_%H-%M")
                pfix = '_psd_fit_results'
                fname = pfix + dstr

        if not fname.endswith(fext):
            fname += fext

        pfile = join(directory, fname)

        cfg = ConfigParser()

        if isfile(pfile):
            cfg.read(pfile)

        for name in self.names:
            fit_results = OrderedDict()

            fit_results[co.model] = self.fits[name].model

            fit_results[co.f_c] = '{0:1.5e}'.format(self.fits[name].f_c)
            fit_results[co.df_c] = '{0:1.5e}'.format(self.fits[name].f_c_err)
            fit_results[co.f_c_unit] = self.fits[name].freq_unit

            fit_results[co.D] = '{0:1.5e}'.format(self.fits[name].D)
            fit_results[co.dD] = '{0:1.5e}'.format(self.fits[name].D_err)
            fit_results[co.D_unit] = self.fits[name].D_unit

            fit_results[co.bounds] = ', '.join(str(a) for a in self.fits[name].get_freq_bounds())
            fit_results[co.redchi2] = '{0:1.5e}'.format(self.fits[name].redchi2)
            fit_results[co.nfree] = self.fits[name].nfree
            fit_results[co.outlier] = self.fits[name].is_outlier()
            fit_results[co.conf_level] = self.fits[name].conf_level

            section = 'FIT_RESULTS_' + name.upper()
            if not section in cfg.sections():
                cfg.add_section(section)
            cfg[section] = fit_results

            # passive calibration results
            pc_results = OrderedDict()
            keys = ['dissens', 'dissens_err', 'dissens_unit',
                    'trap_stiffness', 'trap_stiffness_err',
                    'trap_stiffness_unit',
                    'drag', 'drag_err', 'drag_unit']
            for key in keys:
                key_name = getattr(co, key)
                pc_results[key_name] = getattr(self.pc_results[name], key)

            section = 'PASSIVE_CALIBRATION_RESULTS_' + name.upper()
            if not section in cfg.sections():
                cfg.add_section(section)
            cfg[section] = pc_results

            # active calibration results
            if self.active_calibration:
                ac_results = OrderedDict()

                for key in keys:
                    ac_results[key] = getattr(self.ac_results[name], key)
                section = 'ACTIVE_CALIBRATION_RESULTS_' + name.upper()
                if not section in cfg.sections():
                    cfg.add_section(section)
                cfg[section] = ac_results
        try:
            with open(pfile, 'w') as cfl:
                cfg.write(cfl)
        except:
            raise Exception('Unable to write results to file {}'
                            ''.format(pfile))
        print(pfile)
