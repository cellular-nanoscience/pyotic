# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:50:29 2014

@author: Tobias Jachowski
"""
__author__ = "Tobias Jachowski"
__copyright__ = "Copyright 2016, The PyOTIC Project"
__credits__ = []
__license__ = "Apache-2.0"
__maintainer__ = "Tobias Jachowski"
__email__ = "pyoti@jachowski.de"
__status__ = "beta"


import os
directory = os.path.dirname(globals()['__file__'])

try:
    with open(os.path.join(directory, 'VERSION.txt')) as f:
        __version__ = f.read().strip()
except:
    with open(os.path.join(directory, '..', 'VERSION.txt')) as f:
        __version__ = f.read().strip()


def version():
    """
    Returns the version.
    """
    return __version__


import sys
if 'IPython' in sys.modules:
    try:
        # Check if notebook runs in JupyterLab environment and not the
        # classical notebook
        # TODO: needs to be fixed!
        if False:
            # Try to use ipympl backend (module://ipympl.backend_nbagg)
            import ipympl
        else:
            raise ImportError('Notebook does not run in JupyterLab.')
    except ImportError:
        try:
            ip = get_ipython()

            # Import matplotlib and set backend:
            import matplotlib
            matplotlib.use('nbAgg')
        except:
            pass
            # ip.enable
            # default to inline in kernel environments
            # if hasattr(ip, 'kernel'):
            #     print('enabling inline matplotlib')
            #     ip.enable_matplotlib('inline')
            # else:
            #     print('enabling matplotlib')
            #     ip.enable_matplotlib()

    # Set format for inline plots
    from IPython.display import set_matplotlib_formats
    # %config InlineBackend.figure_formats = ['png']
    set_matplotlib_formats('png', 'svg', 'pdf', 'jpeg', quality=90)

# Load pyoti plugins
from .plugins import plugin_loader
plugin_loader.load_modules()


from . import experiment as ep
from . import evaluate as ev


def info():
    print("PyOTI - the investigator package of the PyOTIC software")
    print("Version: %s" % version())
    print()
    nb_path = os.path.abspath(".")
    print("The actual working path is: '%s'" % nb_path)


def create_experiment(**kwargs):
    """
    Create a new Experiment or get the last created one.
    Optionally, open an experiment file `filename'.
    For documentation see experiment.create_experiment()
    """
    return ep.create_experiment(**kwargs)


def open_experiment(filename=None, directory=None, **kwargs):
    """
    Create a new Experiment or get the last created one.
    Optionally, open an experiment file `filename'.
    For documentation see experiment.open_experiment()
    """
    return ep.open_experiment(filename=filename, directory=directory, **kwargs)


def close_experiment():
    """
    Close the last created Experiment.
    For documentation see experiment.close_experiment()
    """
    return ep.close_experiment()


def save_experiment():
    """
    Save the last created Experiment.
    For documentation see experiment.save_experiment()
    """
    return ep.save_experiment()


def create_calibration(**kwargs):
    """
    Create a Calibration object.
    For documentation see experiment.create_calibration()
    """
    return ep.create_calibration(**kwargs)


def create_tether(**kwargs):
    """
    Create a Tether object.
    """
    return ev.Tether(**kwargs)


def create_motion(**kwargs):
    """
    Create a Motion object.
    """
    return ev.Motion(**kwargs)


def create_stepped(**kwargs):
    raise DeprecationWarning("The function `create_stepped()` was replaced "
                             "by the function `create_motion()`.")
