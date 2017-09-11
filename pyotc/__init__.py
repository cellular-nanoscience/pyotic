# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
__author__ = "Steve Simmert"
__copyright__ = "Copyright 2016, The PyOTIC Project"
__credits__ = []
__license__ = "Apache-2.0"
__maintainer__ = "Steve Simmert"
__email__ = "steve.simmert@uni-tuebingen.de"
__status__ = "stable"

import sys
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


from IPython import get_ipython
import ipykernel

#Load matplotlib and set backend:
ip = get_ipython()
if (hasattr(ipykernel, 'zmqshell')
   and isinstance(ip, ipykernel.zmqshell.ZMQInteractiveShell)):
    import matplotlib
    matplotlib.use('nbAgg')

    from IPython.display import set_matplotlib_formats
    # %config InlineBackend.figure_formats = ['png']
    set_matplotlib_formats('png', 'svg', 'pdf', 'jpeg', quality=90)


import matplotlib.pyplot as plt

# load pint and create unit registry
from pint import UnitRegistry
ureg = UnitRegistry()

from . import focal_shift

from .height_calibration import gen_height_fit_pars
from .height_calibration import fit_rel_drag
from .height_calibration import fit_height_data
from .height_calibration import HeightCalibration
from .height_calibration import HeightFitResult
from .height_calibration import HeightCalibTime

from . import name_constants as co

from .physics import drag
from .physics import faxen_factor

from . import plotting
from .plotting import add_plot_to_figure

from . import psd
from .psd import calculate_psd
from .psd import ExpSetting
from .psd import gen_PSD_from_time_series
from .psd import PSD
from .psd import PSDMeasurement
from .psd import gen_psdm_from_region

from . import psd_fitting
from .psd_fitting import gen_psd_fit_pars
from .psd_fitting import PSDFit

from . import utilities

from scipy.constants import Boltzmann as k_B
