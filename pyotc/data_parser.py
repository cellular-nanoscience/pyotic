# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
"""
Manage the reading and storing of data.
"""

from collections import OrderedDict

from scipy import genfromtxt
from scipy import array

import configparser

from . import name_constants as co


def read_std_data_file(path, lower_names=False, **kwargs):
    '''
    Read a standard data file.

    The file is expected to have the names of the data in the first row. The
    data should be delimited by '\t' (=TABS).

    The output is a dictionary with the names as keys and their respective
    columns as values.
    '''
    data = genfromtxt(path,
                      delimiter='\t',
                      names=True,
                      **kwargs
                      )
    if lower_names:
        names = (name.lower() for name in data.dtype.names)
    else:
        names = data.dtype.names
    cols = array([a for a in zip(*list(data))])

    d = OrderedDict()
    for i, name in enumerate(names):
        d.update({name: cols[i]})

    return d


def read_PSD_parameter_file(path):
    """
    Read a PSD parameter file.

    The parameter file is a chaperone to the PSD data
    file. The file is expected to be in the config file standard.

    The parameters are expected to be in the default section.
    """
    pars = configparser.ConfigParser()

    if not pars.read(path):
        raise Exception('Parameter file could not be read at: {}.'
                        ''.format(path))

    defs = pars.defaults()

    if 'ACTIVE_CALIBRATION' in pars:
        defs.update(pars['ACTIVE_CALIBRATION'])
        defs['active_calibration'] = True
    else:
        defs['active_calibration'] = False

    strings = [co.material, co.ex_axis, 'active_calibration',
               co.names, co.N_avg, co.f_sample, co.medium]

    for k in list(defs):
        if k not in strings and not(k.endswith('unit')):
            defs[k] = float(defs[k])

    return defs


def save_psd_data(path, freq, psd_dict):
    """
    Write psd data to a tab-separated data file.

    Arguments
    ---------
    path : str
        Path to the file.
    freq : array-like
        Frequency vector.
    psd_dict : dict
        Dictionary that holds the data in arrays for each frequency value. The
        Keys are the names of the axes.
    """
    data = OrderedDict()
    data['freq'] = freq
    data.update(psd_dict)
    names = [name for name in data.keys()]

    vals = array([a for a in data.values()]).transpose()

    with open(path, 'w') as fl:
        fl.write('\t'.join(names) + '\n')

        for row in range(len(vals)):
            s = ['{0:1.5E}'.format(si) for si in vals[row]]
            fl.write('\t'.join(s) + '\n')


def save_psd_params(path, param_dict, ac_param_dict=None):
    """
    Saves parameters to a configuration file.

    All key-value pairs in param_dict are written to the 'DEFAULT' section.
    If ac_param_dict dictionary is provided, the key-value pairs are stored
    to the section 'active_calibration'.
    """
    pars = configparser.ConfigParser()
    pars['DEFAULT'] = param_dict

    if ac_param_dict is not None:
        pars['active_calibration'] = ac_param_dict

    with open(path, 'w') as fl:
        pars.write(fl)
