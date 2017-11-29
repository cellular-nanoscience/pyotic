# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 23:31:01 2016

@author: Tobias Jachowski
"""

import configparser
import importlib
import os
import appdirs


pyoti_config_dirs = [os.path.join(os.path.dirname(globals()['__file__']),
                                  'configs')]
site_config_dirs = appdirs.AppDirs('pyoti',
                                   multipath=True).site_config_dir.split(':')
user_config_dirs = appdirs.AppDirs('pyoti',
                                   multipath=True).user_config_dir.split(':')


def print_dirs():
    for dirs in [user_config_dirs, site_config_dirs, pyoti_config_dirs]:
        for d in dirs:
            print(d)


def _cfg_file_paths(cfgfile):
    cfgfiles = []
    for directory in pyoti_config_dirs + site_config_dirs + user_config_dirs:
        file_path = os.path.join(directory, cfgfile)
        if file_path not in cfgfiles:
            cfgfiles.append(file_path)
    return cfgfiles


def read_cfg_file(cfgfile, verbose=True):
    """
    Read in a configuration file and return a dictionary.

    Parameters
    ----------
    cfgfile : str
        The path to the configfile.
    verbose : bool, optional
        Print message, if config file was not found.

    Returns
    -------
    configparser.ConfigParser()
        Includes all sections (keys1) and options (keys2) from the configfile:
        { key1: { key2 : value } }
    """
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfgfiles = _cfg_file_paths(cfgfile)
    cfgfiles_found = cfg.read(cfgfiles)
    if len(cfgfiles_found) == 0:
        raise FileNotFoundError('No such file or directory: %s' % cfgfiles)
    return cfg


def get_cfg_option(cfg, sec, opt, verbose=False):
    """
    Retrieve value of a specific option of a configuration.

    Parameters
    ----------
    cfg : configparser.ConfigParser()
        Configuration as retrieved by the function read_cfg_file().
    sec : str
        The section in which the option is located.
    opt : str
        The option that should be retrieved.
    verbose : bool, optional
        Print info, if either section or option could not be found in cfg.

    Returns
    -------
    str
        Value of the option
    """
    if sec not in cfg:
        verbose and print("Section '%s' is not in configuration '%s'"
                          % (sec, cfg))
        return None

    if opt not in cfg[sec]:
        verbose and print("Option '%s' is not in section '%s'"
                          % (opt, sec))
        return None

    option = cfg[sec][opt]

    return option


def get_cfg_sec_dict(cfg, sec, convert='float', verbose=False):
    """
    Retrieve a dictionary of a section with options as keys and corresponding
    values.

    Parameters
    ----------
    cfg : configparser.ConfigParser()
        Configuration as retrieved by the function read_cfg_file().
    sec : str
        The section in which the options are located.
    convert : str
        The type of the values of the options.
    verbose : bool, optional
        no function so far ...

    Returns
    -------
    dict
        Includes all options (keys) and values from the section selected
        { key: value }
    """
    sec_dic = {}
    if sec in cfg:
        options = cfg[sec]
        for key, value in options.items():
            try:
                if convert == 'int':
                    sec_dic[key] = options.getint(key)
                elif convert == 'float':
                    sec_dic[key] = options.getfloat(key)
                elif convert == 'boolean':
                    sec_dic[key] = options.getboolean(key)
                else:
                    sec_dic[key] = value
            except:
                sec_dic[key] = value

    return sec_dic


def get_cfg_list(cfg, sec, opt, verbose=False):
    """
    Retrieve a comma separated list of values of a specified option.

    Parameters
    ----------
    cfg : configparser.ConfigParser()
        Configuration as retrieved by the function read_cfg_file().
    sec : str
        The section in which the option is located.
    opt : str
        The option that should be retrieved.
    verbose : bool, optional
        Print info, if either section or option could not be found in cfg.

    Returns
    -------
    list
        A list of values of the comma separated option.
    """
    list_ = get_cfg_option(cfg, sec, opt, verbose=verbose)
    if not list_:
        return []
    else:
        return [l.strip() for l in list_.split(',')]


def get_cfg_class(cfg, sec, mod_opt='module', cls_opt='class', std_mod=None,
                  std_cls=None, verbose=False):
    """
    Retrieve a class object as described by an option of a section.

    Parameters
    ----------
    cfg : configparser.ConfigParser()
        Configuration as retrieved by the function read_cfg_file().
    sec : str
        The section in which the option is located.
    mod_opt : str, optional
        The name of the option, which defines the module to load the class from
        (default: 'module').
    cls_opt : str, optional
        The name of the option, wich defines the class to load.
    std_mod : str, optional
        The module to load the class from. Is only used, if mod_opt cannot be
        found in cfg.
    std_cls : str, optional
        The name of the class, if class can not be found
    verbose : bool, optional
        Print info, if section, option or module could not be found.

    Returns
    -------
    class
        The requested class object.
    """
    module_ = get_cfg_option(cfg, sec, mod_opt, verbose=verbose) or std_mod
    if not module_:
        verbose and print("Could not determine module from cfg and no "
                          "std_module given.")
        return None
    try:
        module = importlib.import_module(module_, __package__)
    except:
        verbose and print("Could not load module '%s'" % module_)
        return None

    class_ = get_cfg_option(cfg, sec, cls_opt) or std_cls
    if not class_:
        verbose and print("Could not determine class from cfg.")
        return None

    if not hasattr(module, class_):
        verbose and print("Class '%s' is not in '%s'" % (class_, module))
        return None

    cls = getattr(module, class_)
    return cls
