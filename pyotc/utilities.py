# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """

from . import ureg

from collections import Iterable

from inspect import getargspec

from lmfit import Parameters

from scipy import log10
from scipy import logspace


def gen_fit_pars(**kwargs):
    """
    Generate fit parameters.

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
    params = Parameters()

    for k, v in kwargs.items():
        if k.endswith('_min') or k.endswith('_max'):
            continue

        try:
            v_min = kwargs[k + '_min']
        except:
            v_min = None
        try:
            v_max = kwargs[k + '_max']
        except:
            v_max = None

        params.add(k, value=v, min=v_min, max=v_max)
    return params


def flatten_list(l):
    """
    Flatten an irregular list of lists.

    Found here:
    ('http://stackoverflow.com/questions/2158395/
      flatten-an-irregular-list-of-lists-in-python')
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            for sub in flatten_list(el):
                yield sub
        else:
            yield el


def str2u(s):
    """
    Convert a unit represented as a string into a regular string expression.
    """
    return str(ureg(s).units)


def u2str(unit, formatter='P'):
    """ format a pint unit to a pretty print string"""
    if formatter.lower() == 'p':
        return ''.join('{0:P~}'.format(ureg(unit)).split(' ')[1:])
    elif formatter.lower() == 'l':
        return ''.join('{0:L~}'.format(ureg(unit)).split(' ')[1:])
    else:
        raise ValueError(formatter)


def check_and_separate(fun, kws, exclude=True):
    """
    Checks if keywords of a given dictionary match the signature of
    the given function and separates the keywords.
    """
    argspec = getargspec(fun)

    sep = {}
    for k in kws:
        if k in argspec.args:
            sep[k] = kws[k]

    if exclude:
        for k in sep:
            kws.pop(k)

    return sep


def logspace_points_per_decade(start, end, ppd=5):
    """
    Returns a vector with a fixed number of points per decade.

    When a provided with a range of full decades, like [1, 100] the resulting
    numbers will only differ by a decimal shift. The end will be included.

    Arguments
    ---------
    start : float
        Stating number, e.g. 1.0
    end : float
        last number of the vector, e.g. 1e5
    ppd : int
        Number of points per decade.
    """
    ndec = log10(end) - log10(start)
    return logspace(log10(start), log10(end), num=(ndec * ppd - (ndec-1)))
