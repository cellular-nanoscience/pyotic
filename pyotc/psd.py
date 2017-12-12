# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
"""
Manage power spectral densities measurements.
"""
import pdb

from . import ureg
from . import name_constants as co

from .data_parser import read_PSD_parameter_file
from .data_parser import read_std_data_file
from .data_parser import save_psd_data
from .data_parser import save_psd_params

from .physics import density_H2O
from .physics import drag
from .physics import faxen_factor
from .physics import MATERIALS
from .physics import Material
from .physics import oseen_factor
from .physics import viscosity_H2O
from .physics import dviscosity_H2O

from .plotting import add_plot_to_figure
from .plotting import col_dict

from .utilities import str2u

from collections import OrderedDict

import copy

from inspect import signature

import matplotlib.pyplot as plt

from os.path import join

from scipy import absolute
from scipy import array
from scipy import inf
from scipy import logical_or
from scipy import pi
from scipy import rand
from scipy import randn
from scipy import shape
from scipy import signal
from scipy import sqrt
from scipy import zeros

from scipy.fftpack import fft, ifft, fftfreq

import time

import warnings

#####--------------------------------------------------------------------------
#---- general  ---
#####--------------------------------------------------------------------------

# generate filtered random data
def gen_filtered_data(fun, fs, T_msr, *args, **kwargs):
    """
    Generate filtered random data.

    The function simulates data by calling a random number
    from a normal distribution, which simulates a trapped
    particle with pure white noise. The signal is Fourier-
    transformed and multiplied with the square-root of the given
    filter function. Thus, the filter function must represent the shape
    of a power spectral density. After that it's
    inverse-fourier-transformed.
    The number of data points is fs * T_msr, i.e. the sampling
    frequency times the total measurment time.

    To also account for the pyhsics, one can provide the following two
    keyword arguments:
    mean : float
        Set the mean of the output data.
    std : float
        Set the standard deviation of the output data.

    All other *args and **kwargs are passed to fun(*args, **kwargs).

    Arguments
    ---------
    fun : function
       Filter function, e.g. a low-pass filter.
    fs : float
        Sampling frequency
    T_msr : float
        Measurement time, the total time of the measurement.
    """
    mean = kwargs.pop('mean') if 'mean' in kwargs else 0.0
    std = kwargs.pop('std') if 'std' in kwargs else 1.0

    N = int(fs * T_msr)

    freq = fftfreq(N, 1/fs)

    xraw = randn(N)
    x = ifft(fft(xraw) * sqrt(fun(freq, *args, **kwargs))).real
    if std != 1.0:
        x_out = (x / x.std() * std)
    else:
        x_out = x

    x_out = x_out - x_out.mean() + mean

    return x_out


def calculate_psd(x, fs, N_win=1):
    """
    Calculate power spectral density of the signal x.

    N_win splits the signal into the corresponding number of parts.

    Arguments
    ---------
    x : array(float)
        Signal vector.
    fs : float
        Sampling frequency
    N_win : int
        Number of windows to devide the signal vector into. The power spectra
        are then averaged this many times.
    raise_exception : bool
        Raises an exception instead of a warning if N_win is no common divisor
        of the length of the signal vector.
    """
    N = len(x)
    rest = N % N_win
    if rest > 0:
        warnings.warn('N_win is no common divisor of N = len(x). '
                      '--> {0:1d} data points were omitted.'.format(rest))
        x = x[:-rest]  # throw away some data

    len_ = int(N / N_win)
    psds = []
    for idx in range(N_win):
        freq, psd = signal.welch(x[idx * len_: (idx + 1) * len_],
                                 fs=fs,
                                 window='boxcar',
                                 nperseg=len_,
                                 noverlap=0,
                                 nfft=None,
                                 detrend='constant',
                                 return_onesided=True,
                                 scaling='density',
                                 axis=-1)
        psds.append(psd)

    return (freq, array(psds))


def gen_PSD_from_time_series(x, fs, N_win, calc_errors=False, **PSD_kwargs):
    """
    Generate a **psd** object from a time series data set.

    The function uses the Welch algorithm form the scipy.signal package.

    Arguments
    ---------
    x : array(float)
        Signal vector.
    fs : float
        Sampling frequency.
    N_win : int
        Number of windows to devide the signal vector into. The power spectra
        are then averaged this many times.
    calc_errors : bool
        If True, the errors of the psd values are determined from the single
        (windowed) psds. It is recommended to use the theoretical errors for
        PSD analysis (see Nørrelykke 2010), thus using the default: False.
    PSD_kwargs : keyword arguments
        passed over to the init call of **PSD**.

    Returns
    -------
    PSD
        Power spectral density as an object of the **PSD** class.
    """
    freq, psds = calculate_psd(x, fs, N_win=N_win)

    psd_avg = psds.mean(axis=0)
    if calc_errors:
        err = psds.std(axis=0) / sqrt(N_win)
    else:
        err = psd_avg / sqrt(N_win)

    p = PSD(freq,
            psd_avg,
            err=err,
            N_avg=N_win,
            f_sample=fs,
            **PSD_kwargs)

    return p


#####--------------------------------------------------------------------------
#---- objects ---
#####--------------------------------------------------------------------------


class PSD(object):
    def __init__(self,
                 freq,
                 psd,
                 err=None,
                 name='',
                 f_sample=None,
                 N_avg=1,
                 direction=None,
                 freq_unit='Hz',
                 psd_unit='V**2/Hz'
                 ):
        """
        Describes the one-sided power spectral density for one dimension.

        Arguments
        ---------
        freq : array(float)
            frequency vector of positive frequencies f in **Hertz**. This
            vector could still hold f=0. The psd-value psd(f=0) would then be
            taken as offset.
        psd : array(float)
            Vector with the corresponding power spectral densities. Normally,
            the qunatities are given in Volt²/Hz, since they are, e.g measured
            with a positional sensitive device. Thus the object assumes these
            units.
        err : array(float)
            Corresponding errors to the psd values with the same units.
        name : str
            Gives the PSD a name that is used to categorize it with other PSDs.
            This should idealy be one of 'x' 'y' or 'z' to ensure compatibility
            with HeightCalibration methods, but could be different for other
            purposes, e.g. collective_fit method of PSDFit object, where each
            psd must have a different name.
        f_sample : float
            Sampling frequency.
        N_avg : int
            Number of averages of the psd. In case no error vector is given
            this is  used to  calculate *err*.
        direction : str
            sould be either 'lateral' or 'axial' determining the direction
            in which the fluctuations are measured.
            'axial' means perpendicular to a sample surface.
        freq_unit : str
            Unit of the frequency values ('Hz').
        psd_unit : str
            Unit fo the psd values ('V**2/Hz').
        """
        self.N_avg = N_avg

        self._freq_unit = str2u(freq_unit)
        self._psd_unit = str2u(psd_unit)

        if 0 in freq:
            self.offset = psd[freq == 0][0]
        else:
            self.offset = 0.0

        self._freq = freq[freq > 0]
        self._psd = psd[freq > 0]

        if err is None:
            self._err = self._psd / sqrt(self.N_avg)
        else:
            self._err = err[freq > 0]

        self.reset_mask()

        self.name = name
        if direction is None:
            if 'z' in name:
                self.direction = 'axial'
            else:
                self.direction = 'lateral'
        else:
            self.direction = direction

        if f_sample is None:
            f_sample = 2 * max(freq)
        else:
            if 2 * max(freq) != f_sample:
                warnings.warn('Maximum frequency of the spectrum should, in '
                              'general, be half of the sampling frequency.\n'
                              '2 x f_max = {0:1.3f} != f_sample = {1:1.3f}.\n'
                              'You are able to continue, though.'
                              ''.format(2*max(freq), f_sample))
        self._f_sample = f_sample

    def _exclude_values(self, values, attr_name='freq'):
        """ Excludes values for attr == value """
        attr = getattr(self, attr_name)
        try:
            vals = iter(values)
        except:
            vals = [values]
        for val in vals:
            if val in attr:
                mask = attr == val
                self.mask = logical_or(self.mask, mask)
            else:
                warnings.warn('Value {} not available in attribute '
                              '{}.'.format(val, attr_name))

    def _exclude_values_outside(self, vmin, vmax, attr_name='freq'):
        """
        Exclude values outside the range vmin, vmax of attribute with name
        attr_name'.
        """
        attr = getattr(self, attr_name)
        mask = logical_or(attr < vmin, attr > vmax)
        self.mask = logical_or(self.mask, mask)

    @property
    def freq(self):
        return self._freq[~self.mask]

    def get_freq(self,
                 unit=None,
                 get_all=False,
                 get_masked=False,
                 offset=False):
        """
        Return the frequency vector.

        Arguments
        ---------
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        offset : bool
            Whether to include zero hertz to the freq-vector.
        """
        if unit is None or unit == self._freq_unit:
            conv = 1.0
        else:
            conv = ureg(self._freq_unit).to(unit).magnitude
        if get_masked:
            mask = self.mask
        else:
            mask = ~self.mask
        if get_all:
            f_out = self._freq * conv
        else:
            f_out = self._freq[mask] * conv
        if offset:
            a = list(f_out)
            a.insert(0, 0.0)
            f_out_off = array(a)
            return f_out_off
        else:
            return f_out

    @property
    def psd(self):
        return self._psd[~self.mask]

    def get_psd(self,
                unit=None,
                get_all=False,
                get_masked=False,
                offset=False):
        """
        Return the psd vector.

        Arguments
        ---------
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        offset : bool
            Whether to include the value at zero hertz.
        """
        if unit is None or unit == self._psd_unit:
            conv = 1.0
        else:
            conv = ureg(self._psd_unit).to(unit).magnitude
        if get_masked:
            mask = self.mask
        else:
            mask = ~self.mask
        if get_all:
            p_out = self._psd * conv
        else:
            p_out = self._psd[mask] * conv
        if offset:
            a = list(p_out)
            a.insert(0, self.offset * conv)
            p_out_off = array(a)
            return p_out_off
        else:
            return p_out

    @property
    def psd_err(self):
        return self._err[~self.mask]

    def get_err(self,
                unit=None,
                get_all=False,
                get_masked=False,
                offset=False):
        """
        Return the error vector.

        Arguments
        ---------
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        offset : bool
            Whether to include the error value at zero hertz.
        """
        if unit is None or unit == self._psd_unit:
            conv = 1.0
        else:
            conv = ureg(self._psd_unit).to(unit).magnitude
        if get_masked:
            mask = self.mask
        else:
            mask = ~self.mask
        if get_all:
            e_out = self._err * conv
        else:
            e_out = self._err[mask] * conv
        if offset:
            offset_err = self.offset * conv / sqrt(self.N_avg)
            a = list(e_out)
            a.insert(0, offset_err)
            e_out_off = array(a)
            return e_out_off
        else:
            return e_out

    @property
    def f_sample(self):
        return self._f_sample

    def get_f_sample(self, unit=None):
        if unit is None or unit == self._freq_unit:
            return self._f_sample
        else:
            conv = ureg(self._freq_unit).to(unit).magnitude
            return self._f_sample * conv

    def set_f_sample(self, f_sample, unit):
        if unit == self._f_sample:
            conv = 1.0
        else:
            conv = ureg(unit).to(self._freq_unit).magnitude
        self._f_sample = f_sample * conv

    @property
    def df(self):
        """ Frequency resolution of the power spectrum"""
        return self._freq[1] - self._freq[0]

    @property
    def T_msr(self):
        """ Measurement time in seconds """
        return (1 / min(self.get_freq(unit='Hz', get_all=True)))

    @property
    def N_samples(self):
        """ Number of samples """
        return (self.T_msr * self.f_sample)

    def is_lateral(self):
        """
        Whether the direction attribute is set to 'lateral' (True) or
        'axial' (False).
        """
        if self.direction == 'lateral':
            return True
        elif self.direction == 'axial':
            return False
        else:
            raise Exception('Unknown direction {}'.format(self.direction))

    def reset_mask(self):
        """ Rest the internal mask, so all values a taken into account. """
        self.mask = zeros(shape(self._freq)) > 0

    def add_mask(self, mask):
        """
        Logical OR with self.mask.
        """
        self.mask = logical_or(self.mask, mask)

    def exclude_freq(self, f_exclude):
        """
        Exclude the data at frequency or frequencies. f_exclude can be both a
        list of or a single float.
        """
        self._exclude_values(f_exclude, attr_name='_freq')

    def exclude_freq_outside(self, fmin, fmax):
        """
        Exclude the frequencies outside the interval [fmin, fmax].
        """
        self._exclude_values_outside(fmin, fmax, attr_name='_freq')

    def plot_psd(self,
                 axis=None,
                 plot_all=False,
                 plot_masked=False,
                 plot_errors=False,
                 **kwargs
                 ):
        """
        Plots the power spectral density.

        Arguments
        ---------
        axis : axis
            Axis to plot to.
        plot_all : bool
            If True, omits the internal mask.
        plot_masked : bool
            When True plots the masked values instead.
        plot_errors : bool
            Plot error bars as well.
        **kwargs
            key-word arguments handed over to plot()

        Returns
        -------
        figure
        """
        freq = self.get_freq(get_all=plot_all, get_masked=plot_masked)
        psd = self.get_psd(get_all=plot_all, get_masked=plot_masked)
        err = self.get_err(get_all=plot_all, get_masked=plot_masked)

        if plot_masked:
            if self.name in col_dict:
                col = col_dict['o' + self.name]
            else:
                col = 'gray'
            fmt = 'o'
            if 'markersize' not in kwargs:
                kwargs['markersize'] = 3
        else:
            if self.name in col_dict:
                col = col_dict[self.name]
            else:
                col = tuple(rand(3))
            fmt = '-'

        if not plot_errors:
            err = []

        if axis is None:
            fig = None
        else:
            fig = axis.figure

        if 'color' not in kwargs and self.name in col_dict.keys():
            kwargs['color'] = col
        if 'fmt' not in kwargs:
            kwargs['fmt'] = fmt
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 1.0
        if 'alpha' not in kwargs:
            kwargs['alpha'] = 0.5
        if 'title' not in kwargs:
            kwargs['title'] = 'PSD {}'.format(self.name)
        if 'showLegend' not in kwargs:
            kwargs['showLegend'] = True
        if 'legend_kwargs' not in kwargs:
            lg_kws = {'loc': 3}
            kwargs['legend_kwargs'] = lg_kws
        if 'fontsize' not in kwargs:
            kwargs['fontsize'] = 16

        ax = add_plot_to_figure(fig,
                                freq,
                                psd,
                                yerr=err,
                                label=self.name,
                                axis=axis,
                                logplot=True,
                                **kwargs
                                )
        ax.grid(which='major')
        plt.setp(ax,
                 xlabel='Frequency (Hz)',
                 ylabel=r'PSD ($\mathsf{V^2/Hz}$)'
                 )

        return ax.figure


class ExpSetting(object):
    """
    Describes the experimental setting of a Measurment.
    """
    def __init__(self,
                 temp,
                 radius,
                 temp_err=0.0,
                 radius_err=0.0,
                 height=inf,
                 density_particle=None,
                 density_medium=None,
                 viscosity=None,
                 viscosity_err=0.0,
                 material='',
                 medium='water',
                 temp_unit='K',
                 radius_unit='m',
                 height_unit='m',
                 density_particle_unit='kg/m**3',
                 density_medium_unit='kg/m**3',
                 viscosity_unit='Pa*s',
                 warn=True):
        """
        Define the experimental setting.

        Arguments
        ---------
        temp : float
            Temperature at which the measurement was performed.
        radius : float
            Radius of the bead.
        temp_err : float
            Error of the Temperature value.
        radius_err : float
            Radius error.
        height : float
            Apperent distance between trapped particle and surface.
        density_particle : float
            The mass density of the bead. If this is provided, **material**
            will only name the material, but has no other effect.
        density_medium : float or callable.
            Density of the medium. If callable, then the function is called
            as density_medium(temp), with temp in Kelvin. Otherwise the float
            value is used.
        viscosity : float or callable
            Viscosity of the medium. If callable, then the function is called
            as viscosity(temp), with temp in Kelvin. Otherwise the float value
            is used.
        viscosity_err : float, callable or None
            If callable, the function is used as
            absolute(viscosity_err(temp) * temp_err), whit temp in Kelvin.
            Otherwise the float value is used as error in the viscosity
            estimate. If it's None, a known medium must be given.
        material : str
            Name of the particle material. E.g. 'polystyrene', 'titania',
            'silica'. If the name is none of these the particle density must be
            provided.
        temp_unit : str
            Unit of the Temperature ('K').
        radius_unit : str
            Unit of the radius values (e.g. 'um').
        height_unit : str
            Unit of the height.
        density_particle_unit : str
            The unit the density is provided in, e.g. 'kg/m**3'
        density_medium_unit : str
            The unit the density is provided in, e.g. 'kg/m**3'
        viscosity_unit : str
            Unit of the provided viscosity, e.g. 'Pa*s'.

        Notes
        -----
        material
            If material is a known material a set desnity_particle will not
            have an effect.
        copy
            Use copy() to make a duplicate of the setting.
        """
        # temperature
        self._qtemp = ureg.Measurement(temp, temp_err, temp_unit)

        # radius
        self._radius = radius
        self._radius_err = radius_err
        self._radius_unit = str2u(radius_unit)

        # height
        self._height = height
        self._height_unit = height_unit

        # material, density_material
        try:
            self.material = Material(material,
                                     density=density_particle,
                                     density_unit=density_particle_unit)
        except:
            if density_particle:
                if material == '':
                    material_name = 'Unknown material'
                else:
                    material_name = material
                self.material = Material(material_name,
                                         density=density_particle,
                                         density_unit=density_particle_unit)
            else:
                if warn:
                    warnings.warn('Unknown material or density not given. '
                                  'Using fallback density=1050 kg/m**3 '
                                  '(Polystyrene).')
                self.material = Material('PS')
        self.get_density_particle = self.material.get_density

        # medium
        if medium.lower() in ['water', 'h2o']:
            self._density_medium = density_H2O(self.get_temp(unit='K'))
            self._density_medium_unit = str2u('kg/m**3')
            self._viscosity = viscosity_H2O(self.get_temp(unit='K'))
            dv = absolute(dviscosity_H2O(self.get_temp(unit='K')) *
                          self.get_temp_err(unit='K'))
            self._viscosity_err = dv
            self._viscosity_unit = str2u('Pa*s')
            self.medium = 'water'
        else:
            if density_medium and viscosity:
                T = self.get_temp(unit='K')
                dT = self.get_temp_err(unit='K')
                try:
                    self._density_medium = density_medium(T)
                except:
                    self._density_medium = density_medium
                self._density_medium_unit = str2u(density_medium_unit)
                try:
                    self._viscosity = viscosity(T)
                except:
                    self._viscosity = viscosity
                try:
                    self._viscosity_err = absolute(viscosity_err(T) * dT)
                except:
                    self._viscosity_err = viscosity_err
                self._viscosity_unit = str2u(viscosity_unit)
                self.medium = medium
            else:
                raise Exception('You need to specify density and viscosity of '
                                'the medium "{}"'.format(medium))

    @property
    def _temp_unit(self):
        return str(self._qtemp.units)

    @property
    def temp(self):
        return self._qtemp.value.magnitude

    def get_temp(self, unit=None):
        if unit is None or unit == self._temp_unit:
            return self.temp
        else:
            return self._qtemp.to(unit).value.magnitude

    def set_temp(self, temp, unit):
        """ Set the Temperature of the measurement in specified units. """
        T_err = self._qtemp.error.magnitude
        if unit == self._temp_unit:
            self._qtemp = ureg.Measurement(temp, T_err, unit)
        else:
            T = ureg.Measurement(temp, T_err, unit)
            self._qtemp = T.to(self._temp_unit)

    @property
    def temp_err(self):
        return self._qtemp.error.magnitude

    def get_temp_err(self, unit=None):
        if unit is None or unit == self._temp_unit:
            return self.temp_err
        else:
            return self._qtemp.to(unit).error.magnitude

    def set_temp_err(self, T_err, unit):
        """ Set the Temperature of the measurement in specified units. """
        T = self.temp
        if unit == self._temp_unit:
            self._qtemp = ureg.Measurement(T, T_err, unit)
        else:
            T_ = ureg.Measurement(T, T_err, unit)
            self._qtemp = T_.to(self._temp_unit)

    @property
    def radius(self):
        return self._radius

    def get_radius(self, unit=None):
        if unit is None or unit == self._radius_unit:
            return self._radius
        else:
            conv = ureg(self._radius_unit).to(unit).magnitude
            return self._radius * conv

    def set_radius(self, radius, unit):
        """
        Set the radius of the used bead in the specified unit.

        Unit is then converted into the preset units.
        """
        if unit == self._radius_unit:
            conv = 1.0
        else:
            conv = ureg(unit).to(self._radius_unit).magnitude
        self._radius = radius * conv

    @property
    def radius_err(self):
        return self._radius_err

    def get_radius_err(self, unit=None):
        if unit is None or unit == self._radius_unit:
            return self._radius_err
        else:
            conv = ureg(self._radius_unit).to(unit).magnitude
            return self._radius_err * conv

    def set_radius_err(self, radius_err, unit):
        """
        Set the error of the radius of the used bead in the specified unit.

        Unit is then converted into the preset units.
        """
        if unit == self._radius_unit:
            conv = 1.0
        else:
            conv = ureg(unit).to(self._radius_unit).magnitude
        self._radius_err = radius_err * conv

    @property
    def height(self):
        """ Apperent distance between trapped particle and surface."""
        return self._height

    def get_height(self, unit=None):
        """
        Return the apperent distance between trapped particle and surface in
        the specified units.
        """
        if unit is None or unit == self._height_unit:
            return self._height
        else:
            conv = ureg(self._height_unit).to(unit).magnitude
            return self._height * conv

    def set_height(self, height, unit):
        """
        Set the apperent distance between trapped particle and surface in
        specified units.
        """
        conv = ureg(unit).to(self._height_unit).magnitude
        self._height = height * conv

    @property
    def density_particle(self):
        """Return the particle density in kg/m**3."""
        return self.material.get_density(unit='kg/m**3')

    #NOTE: get_density_particle is defined in __init__ and is hooked to
    #      material.get_density

    def set_material(self, name, density, density_unit):
        self.material = Material(name,
                                 density=density,
                                 density_unit=density_unit)

    @property
    def viscosity(self):
        return self._viscosity

    def get_viscosity(self, unit=None):
        if unit is None or unit == self._viscosity_unit:
            conv = 1.0
        else:
            conv = ureg(self._viscosity_unit).to(unit).magnitude
        return self._viscosity * conv

    def set_viscosity(self, viscosity, unit='Pa*s'):
        """
        Set the viscosity.

        Arguments
        ---------
        viscosity : float or callable
            Viscosity of the medium. If callable, then the function is called
            as viscosity(temp), with temp in Kelvin. Otherwise the float value
            is used.
        unit : str
            unit of the viscosity
        """
        T = self.get_temp(unit='K')
        try:
            self._viscosity = viscosity(T)
        except:
            self._viscosity = viscosity
        self._viscosity_unit = str2u(unit)

    def get_viscosity_err(self, unit=None):
        if unit is None or unit == self._viscosity_unit:
            conv = 1.0
        else:
            conv = ureg(self._viscosity_unit).to(unit).magnitude
        return self._viscosity_err * conv

    def set_viscosity_err(self, viscosity_err, unit='Pa*s'):
        """
        Set the viscosity.

        Arguments
        ---------
        viscosity_err : float or callable
            If callable, the function is used as
            absolute(viscosity_err(temp) * temp_err), whit temp in Kelvin.
            Otherwise the float value is used as error in the viscosity
            estimate.
        unit : str
            unit of the viscosity
        """
        T = self.get_temp(unit='K')
        dT = self.get_temp_err(unit='K')
        try:
            self._viscosity_err = absolute(viscosity_err(T) * dT)
        except:
            self._viscosity_err = viscosity_err
        self._viscosity_unit = str2u(unit)

    @property
    def density_medium(self):
        return self._density_medium

    def get_density_medium(self, unit=None):
        if unit is None or unit == self._density_medium_unit:
            conv = 1.0
        else:
            conv = ureg(self._density_medium_unit).to(unit).magnitude
        return self._density_medium * conv

    def set_density_medium(self, density, unit='Pa*s'):
        self._density_medium = density
        self._density_medium_unit = str2u(unit)

    def copy(self):
        return copy.copy(self)

    def get_dict(self):
        """
        Return data inside a dictionary.
        """
        d = OrderedDict()

        d[co.temp] = self.temp
        d[co.temp_err] = self.temp_err
        d[co.temp_unit] = self._temp_unit

        d[co.radius] = self.radius
        d[co.radius_err] = self.radius_err
        d[co.radius_unit] = self._radius_unit

        d[co.height] = self.height
        d[co.height_unit] = self._height_unit

        d[co.material] = self.material.name
        if self.material.name not in MATERIALS:
            d[co.density_p] = self.density_particle
            d[co.density_p_unit] = str2u('kg/m**3')

        d[co.medium] = self.medium
        if self.medium != 'water':
            d[co.density_m] = self.density_medium
            d[co.density_m_unit] = self._density_medium_unit

            d[co.viscosity] = self.viscosity
            d[co.viscosity_unit] = self._viscosity_unit

        return d


class PSDMeasurement(object):
    def __init__(self,
                 exp_setting=None,
                 freq_unit='Hz',
                 psd_unit='V**2/Hz',
                 warn=True):
        """
        Manage power spectra (PSDs) of one measurement.

        One-sided PSDs are assumed! They are stored in self.psds.
        Use **load()** to read data from a PSD data file or
        **add_psd()** to add **psd**-objects one after another.

        Arguments
        ---------
        exp_setting : ExpSetting
            Object discribing the experimental setting. Satndard values are
            assumed, if None. A warning will show up if warn is not False.
        freq_unit : str
            Unit of the frequency values ('Hz').
        psd_unit : str
            Unit fo the psd values ('V**2/Hz').
        warn : bool
            Show warnings if standard values are assumend.
        """
        if not exp_setting:
            temp = 25.0
            radius = 0.5
            exp_setting = ExpSetting(temp, radius,
                                     temp_unit='degC',
                                     radius_unit='um',
                                     warn=warn)
            if warn:
                warnings.warn('Standard values for experimental setting used! '
                              'I.e. temp=25 degC, radius=0.5 um, height=inf,'
                              'material=Polystyrene')

        if not isinstance(exp_setting, ExpSetting):
            raise Exception('exp_setting is not an ExpSetting object.')

        self.exp_setting = exp_setting

        self.psds = OrderedDict()

        self._freq_unit = freq_unit
        self._psd_unit = psd_unit

        self.active_calibration = False
        self.ex_axis = None
        self.ex_freq = None
        self.ex_amplitude = None
        self.ex_amplitude_err = None
        self._ex_amplitude_unit = None

        self.ex_power = None
        self.ex_power_err = None

        self._ex_power_unit = None

    def add_psd(self, name, psd):
        """
        Add a **psd** object to the measurement.

        Arguments
        ---------
        name : str
            String specifying the name of the axis.
        psd : PSD
            Object of the **PSD** class
        """
        # The first psd that's added set's the standard values of other psd
        # that might be added later on.
        if name in self.get_names():
            warnings.warn("PSD-object with name '{}' already "
                          "present. Data overridden".format(name))
        psd.name = name
        self.psds[name] = psd

    @property
    def names(self):
        """names of the axes."""
        return self.get_names()

    def get_names(self):
        """Return the names of the axes."""
        lst = [psd.name for psd in self.psds.values()]
        lst.sort()
        return lst

    def get_f_sample(self, name, unit=None):
        """
        Return fampling frequency of the PSD with the given name.
        """
        if unit is None:
            unit = self._freq_unit
        return self.psds[name].get_f_sample(unit=unit)

    def get_f_resolution(self, name, unit=None):
        if unit is None:
            unit = self._freq_unit
        unit_psd = self.psds[name]._freq_unit
        if unit_psd != unit:
            conv = ureg(unit_psd).to(unit).magnitude
        else:
            conv = 1.0
        return self.psds[name].df * conv

    def get_N_samples(self, name):
        return self.psds[name].N_samples

    def get_laterality(self):
        """
        Return a dictionary with names of the psds as keys and values
        discribing if the direction is lateral (True) or axial (False).
        """
        d = {name: psd.is_lateral()
             for name, psd in self.psds.items()}
        return d

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
        Return the frequency vector of the axis specified by **name**.

        Arguments
        ---------
        name : str
            name of the axis.

        Keyword Arguments
        -----------------
        unit : str
        get_all : bool
            If False uses the internal mask.
        get_masked : bool
            If True returns the masked instead of the unmasked elements.
        """
        return self.psds[name].get_freq(**kwargs)

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
        Return the psd vector of the axis specified by 'name'.

        Arguments
        ---------
        name : str
            name of the axis.
        kwargs
            can be:
             unit : str
             get_all : bool
                 If False uses the internal mask.
             get_masked : bool
                 If True returns the masked instead of the unmasked elements.
        """
        return self.psds[name].get_psd(**kwargs)

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
        Return the error vector of the axis specified by 'name'.

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
        offset : bool
            Whether to include the value at zero hertz.
        """
        return self.psds[name].get_err(**kwargs)

    def set_ac_params(self,
                      ex_axis,
                      ex_freq,
                      ex_amplitude,
                      ex_power,
                      ex_amplitude_err=0,
                      ex_power_err=0,
                      freq_unit='Hz',
                      amplitude_unit='m',
                      power_unit='V**2'):
        """
        Set the values of an active PSD measurement.

        Arguments
        ---------
        ex_axis : str
            Defines the name of the psd that was excited this must be one
            of self.names.
        ex_freq : float
            Frequency at which the bead was driven.
        ex_amplitude : float
            Amplitude of the sine movement of the bead w.r.t. the trap center.
        amplitude_unit : str
            Unit that the amplitude is given in. The amplitude is then
            converted to the same unit as the radius unit.
        ex_power : float
            The experimentally determined power at the given frequency.
        power_unit : str
            Unit the power is given in. This should be equal to
            psd unit * freq unit.
        ex_amplitude_err : float
            Error of the amplitude.
        ex_power_err : float
            Error of the measured power.

        Note
        ----
        It is good practise to measure the driving amplitude, e.g. via the
        piezo stage monitor signal, instead of taking the set amplitude for
        granted. The set value often differs from the actual one,
        especially at high driving frequencies.

        Reference
        ---------
        Tolić-Nørrelykke, et al. (2006)
            Part IV A
            Calibration of optical tweezers with positional detection in
            the back focal plane.
            Review of Scientific Instruments, 77(10), 103101.
            http://doi.org/10.1063/1.2356852
        """
        if ex_axis not in self.names:
            raise Exception('Unknown axis {}'.format(ex_axis))
        else:
            self.ex_axis = ex_axis

        if freq_unit != self._freq_unit:
            conv = ureg(freq_unit).to(self._freq_unit).magnitude
        else:
            conv = 1.0
        self.ex_freq = float(ex_freq) * conv

        r_unit = self.exp_setting._radius_unit
        if amplitude_unit != r_unit:
            conv = ureg(amplitude_unit).to(r_unit).magnitude
        else:
            conv = 1.0
        self.ex_amplitude = float(ex_amplitude) * conv
        self.ex_amplitude_err = float(ex_amplitude_err) * conv
        self._ex_amplitude_unit = r_unit

        self.ex_power = float(ex_power)
        self.ex_power_err = float(ex_power_err)
        self._ex_power_unit = power_unit

        self.active_calibration = True

    def get_stokes_drag(self, unit='N*s/m'):
        """
        Return Stoke's drag in the specified unit.

        The function uses the drag function of pyotc.physics evaluated at
        f=0 and height=inf.

        See Also
        --------
        pyotc.physics.drag
        """
        d = drag(self.exp_setting.get_radius(unit='m'),
                 self.exp_setting.get_temp(unit='K'),
                 density=self.exp_setting.get_density_medium(unit='kg/m**3'),
                 viscosity=self.exp_setting.get_viscosity(unit='Pa*s')).real
        if unit == 'N*s/m':
            conv = 1.0
        else:
            conv = ureg('N*s/m').to(unit).magnitude
        return d * conv

    def get_stokes_drag_err(self, unit='N*s/m'):
        """
        Return the error to Stoke's drag in the specified unit.

        The function calculates the error on the calculated stokes drag from
        the relative errors of the viscosity and the radius.

        See Also
        --------
        get_stokes_drag
        """
        v = self.exp_setting.get_viscosity()
        dv = self.exp_setting.get_viscosity_err()

        dr = self.exp_setting.get_radius_err(unit='m')
        r = self.exp_setting.get_radius(unit='m')
        ddrag = (dv / v + dr / r) * self.get_stokes_drag(unit=unit)
        return ddrag

    def get_corrected_drag(self,
                           mode=1,
                           distance=None,
                           focal_shift=1.0,
                           distance_unit='um',
                           drag_unit='N*s/m'):
        """
        Correct stokes drag for presence of either one wall or two parallel
        walls at the given height of the trapped particle.

        Only parallel movement of the sphere is considered!

        mode = 1
            Faxén's correction: drag = drag_stokes * faxen_factor
        mode = 2
            Linear superposition approach as done by Oseen in 1927; you need
            to provide the thickness of the sample!

            drag = drag_stokes * (faxen_factor(h) + faxen_factor(H-h)-1)

        Arguments
        ---------
        mode : 1 or 2
            If one or two walls are present.
        distance : float
            The distance between two parallel walls.
        focal_shift : float defaults to 1.0
            Relative shift of the focal point of the trap, due to diffraction.
            This usually is about 0.8 for an oil-immersion objective and a
            sample in water. The given height is automatically corrected for
            that shift.
        distance_unit : str
            Unit of distance.
        drag_unit : str
            unit of the output drag.
        """
        drag_stokes = self.get_stokes_drag(unit=drag_unit)
        radius = self.exp_setting.get_radius(unit=distance_unit)
        height = self.exp_setting.get_height(unit=distance_unit) * focal_shift

        if height < 1.5 * radius:
            warnings.warn("Faxén's correction does not work for heights lower"
                          "than 1.5 * radius.")

        if mode == 1:
            drag = drag_stokes * faxen_factor(height, radius)
        elif mode == 2:
            if distance is None:
                raise Exception('Distance between the two walls is not given.')
            if height > (distance - 1.5 * radius):
                warnings.warn('The correction does not work for heights above'
                              'thickness - 1.5 * radius.')
            drag = drag_stokes * oseen_factor(height, radius, distance)
        else:
            raise Exception('Unknown mode {}'.format(mode))
        return drag

    def save(self,
             directory=None,
             datafile=None,
             suffix='_psd_parameters.txt',
             datafile_extension='.dat',
             include_errors=False
             ):
        """
        Save PSD measurment to data and parameter file.

        Arguments
        ---------
        directory : path
            If None is given and there is no attribute 'directory' the current
            working directory is used.
        datafile : str
            If None, and there is no attribute 'datafilename' the current date
            and time is used.
        suffix : str
            Suffix for the parameter file.
        datafile_extension : str
            Extension of the data file.
        """
        # save data in data file
        if directory is None:
            try:
                directory = self.directory
            except:
                directory = './'
        if datafile is None:
            try:
                datafile = self.datafilename
            except:
                dstr = time.strftime("%Y-%m-%d_%H-%M")
                pfix = '_psd_measurement_'
                datafile = pfix + dstr

        if not datafile.endswith(datafile_extension):
            datafile += datafile_extension

        ptdfile = join(directory, datafile)

        # check that frequency vectors have the same lengths
        # TODO alternative: save freq vectors with axis name 'freq_x' etc.
        flengths = [len(self.get_freq(ax, get_all=True, offset=True))
                    for ax in self.names]
        if len(set(flengths)) > 1:
            raise Exception('Frequency vectors have different lengths!')

        psd_dict = OrderedDict()
        psd_dict.update({'psd_' + ax: self.get_psd(ax,
                                                   get_all=True,
                                                   offset=True)
                         for ax in self.names})
        if include_errors:
            psd_dict.update({'err_' + ax: self.get_psd_err(ax,
                                                           get_all=True,
                                                           offset=True)
                             for ax in self.names})

        plengths = [len(psd) for psd in psd_dict.values()]
        if len(set(plengths)) > 1:
            raise Exception('PSD vectors have different lengths!')

        freq = self.get_freq(self.names[0], get_all=True, offset=True)

        save_psd_data(ptdfile, freq, psd_dict)

        # save parameters in parameter file
        pfile = datafile[:datafile.rfind(datafile_extension)] + suffix
        ptpfile = join(directory, pfile)

        # generate parameters for parameter file
        params = OrderedDict()
        # names
        params[co.names] = ','.join(self.names)
        # N_avg
        N_avg_ = [self.psds[name].N_avg for name in self.names]
        params[co.N_avg] = ','.join(str(n) for n in N_avg_)
        # f_sample
        fs_ = [self.get_f_sample(name, unit=self._freq_unit)
               for name in self.names]
        params[co.f_sample] = ','.join(str(f) for f in fs_)

        params[co.freq_unit] = self._freq_unit
        params[co.psd_unit] = self._psd_unit

        params.update(self.exp_setting.get_dict())

        if self.active_calibration:
            ac_params = OrderedDict()
            ac_params[co.ex_axis] = self.ex_axis
            ac_params[co.ex_freq] = self.ex_freq
            ac_params[co.ex_amp] = self.ex_amplitude
            ac_params[co.ex_amp_err] = self.ex_amplitude_err
            ac_params[co.ex_pow] = self.ex_power
            ac_params[co.ex_pow_err] = self.ex_power_err
            ac_params[co.ex_amp_unit] = self.exp_setting._radius_unit
            ac_params[co.ex_pow_unit] = self._ex_power_unit
        else:
            ac_params = None

        save_psd_params(ptpfile, params, ac_param_dict=ac_params)

        self.directory = directory
        self.datafilename = datafile
        self.paramfile = pfile

    def load(self,
             directory,
             datafile,
             paramfile=None,
             suffix='_psd_parameters.txt'
             ):
        """
        Load the psd-data and according parrameter file.

        Arguments
        ---------
        directory : path
            Path to the data folder.
        datafile : str
            Filename of the psd-data file '*.dat'.
        paramfile : str
            Filename of the corresponding parameter file. If None, the
            extensiotn of the datafile is cut off and **suffix** is appended
        suffix : str
            Suffix to the datafile name that characterizes the parameter file.

        Note
        ----
        datafile
            The data file is assumed to have the data in columns, with their
            names in the first row. Normally the columns are called 'freq',
            'psd_x', psd_y, etc. Columns starting with 'psd_' are still loaded,
            and the strings after the underscore '_' are taken as names for the
            created **psd** objects. So 'PSD_x1', 'PSD_x2', etc. will also
            work and create axes with names x1, x2, etc.
            **Important:** Any name that contains 'z' is considered to be an
            axis in the axial direction. Hence **psd.direction** is set to
            'axial'.
        paramfile
            The parameter file is read by the module **configparser**,
            so a **config-file standard** is assumed.
            Parameters of the [DEFAULT] section are read.
            The following parameters are read, the values in brackets are the
            fallback values - if none is given the value is mendatory:
             - freq_unit ('Hz')
             - n_avg
             - sampling_rate or f_sample
             - psd_unit ('V**2/Hz')
             - bead_dia or radius
             - err_bead_dia or radius_err (0.0)
             - diameter_unit or radius_unit ('um')
             - density and density_unit or material ('PS')
             - density_med (None - water is used later on)
             - density_med_unit ('kg/m**3')
             - viscosity (None - water is used later on)
             - viscosity_unit ('Pa*s')
             - height (inf)
             - height_unit ('um')
             - temperature (25)
             - temperature_error (5)
             - temp_unit ('celsius')
            If the trapped particle was actively driven by a sine-wave and
            an active calibration should be done, the following values must be
            in the 'ACTIVE_CALIBRATION' section of the parameter file.
             - excitation_axis
             - excitation_frequency
             - excitation_amplitude
             - excitation_amplitude_error (0)
             - amplitude_unit (um)
             - reference power or power
             - reference_power_error or power_err
             - power_unit ('V**2')
        """
        dfile = join(directory, datafile)

        if paramfile is None:
            paramfile = datafile[:-4] + suffix

        self.paramfile = paramfile

        pfile = join(directory, paramfile)

        # read psd data file
        data = read_std_data_file(dfile, lower_names=True)

        # read parameter file
        pars = read_PSD_parameter_file(pfile)

        # set up fallback values
        params = {co.freq_unit: 'Hz',
                  co.psd_unit: 'V**2/Hz',
                  co.radius_err: 0.0,
                  co.radius_unit: 'm',
                  co.height: inf,
                  co.height_unit: 'm',
                  co.temp: 25.0,
                  co.temp_err: 2.5,
                  co.temp_unit: 'celsius',
                  co.material: '',
                  co.density_p: None,
                  co.density_p_unit: 'kg/m**3',
                  co.medium: 'water',
                  co.density_m: None,
                  co.density_m_unit: 'kg/m**3',
                  co.viscosity: None,
                  co.viscosity_unit: 'Pa*s'}
        # overwrite default parameters
        params.update(pars)

        # backward compatibility to version 0.2.2
        if co.names not in params:
            # get the axis names from the header in the datafile
            names = [name.split('_')[1]
                     for name in data.keys() if name.lower().startswith('psd')]
        else:
            names = params[co.names].split(',')

        # check if N_avg and f_sample are specified
        for par in [co.N_avg, co.f_sample]:
            if par not in params:
                raise Exception('Parameter {} missing in parameter file {}'
                                ''.format(par, pfile))

        n_avg_ = params[co.N_avg].split(',')
        # backward compatibility to version 0.2.2
        if len(n_avg_) == 1:
            n_avg_ = [n_avg_[0] for name in names]
        N_avg = {}
        for name, n_avg in zip(names, n_avg_):
            N_avg[name] = int(n_avg)

        f_sample_ = params[co.f_sample].split(',')
        # backward compatibility to version 0.2.2
        if len(f_sample_) == 1:
            f_sample_ = [f_sample_[0] for name in names]
        f_sample = {}
        for name, fs in zip(names, f_sample_):
            f_sample[name] = float(fs)

        # check if radius or diameter was specified
        if co.radius not in params:
            raise Exception('Bead radius is not specified in '
                            'parameter file {}'.format(paramfile))

        expset = ExpSetting(params[co.temp],
                            params[co.radius],
                            temp_err=params[co.temp_err],
                            radius_err=params[co.radius_err],
                            height=params[co.height],
                            density_particle=params[co.density_p],
                            density_medium=params[co.density_m],
                            viscosity=params[co.viscosity],
                            material=params[co.material],
                            medium=params[co.medium],
                            temp_unit=params[co.temp_unit],
                            radius_unit=params[co.radius_unit],
                            height_unit=params[co.height_unit],
                            density_particle_unit=params[co.density_p_unit],
                            density_medium_unit=params[co.density_m_unit],
                            viscosity_unit=params[co.viscosity_unit])
        self.exp_setting = expset

        # now get teh freq vector and psd data and put it into a PSD object
        freq = data.pop(co.freq)

        for name in names:
            try:
                psd_vals = data['psd_' + name]
            except:
                warnings.warn('PSD value for axis name {} not found'
                              ''.format(name))
                continue
            try:
                psd_err = data['err_' + name]
            except:
                psd_err = None

            psd = PSD(freq,
                      psd_vals,
                      err=psd_err,
                      name=name,
                      f_sample=f_sample[name],
                      N_avg=N_avg[name],
                      freq_unit=params[co.freq_unit],
                      psd_unit=params[co.psd_unit]
                      )
            self.add_psd(name, psd)

        # check if there's active calibration information
        if params['active_calibration']:
            ex_axis = params[co.ex_axis]
            ex_freq = params[co.ex_freq]
            ex_amplitude = params[co.ex_amp]
            ex_power = params[co.ex_pow]

            if co.ex_amp_err in params:
                ex_amplitude_err = params[co.ex_amp_err]
            else:
                ex_amplitude_err = 0

            if co.ex_pow_err in params:
                ex_power_err = params[co.ex_pow_err]
            else:
                ex_power_err = 0

            if co.ex_amp_unit in params:
                ex_amp_unit = params[co.ex_amp_unit]
            else:
                ex_amp_unit = self._radius_unit

            if co.ex_pow_unit in params:
                power_unit = params[co.ex_pow_unit]
            else:
                psd_unit = params[co.psd_unit]
                freq_unit = params[co.freq_unit]
                power_unit = str((ureg(psd_unit) * ureg(freq_unit)).units)

            self.set_ac_params(ex_axis,
                               ex_freq,
                               ex_amplitude,
                               ex_power,
                               ex_amplitude_err=ex_amplitude_err,
                               ex_power_err=ex_power_err,
                               freq_unit=params[co.freq_unit],
                               amplitude_unit=ex_amp_unit,
                               power_unit=power_unit)

        self.directory = directory
        self.datafilename = datafile

    def get_ex_freq(self, unit='Hz'):
        if unit != self._freq_unit:
            conv = ureg(self._freq_unit).to(unit).magnitude
        else:
            conv = 1.0
        return self.ex_freq * conv

    def get_ex_amplitude(self, unit=None):
        if unit is None or unit == self._ex_amplitude_unit:
            return self.ex_amplitude
        else:
            conv = ureg(self._ex_amplitude_unit).to(unit).magnitude
            return self.ex_amplitude * conv

    def get_ex_amplitude_err(self, unit=None):
        if unit is None or unit == self._ex_amplitude_unit:
            return self.ex_amplitude_err
        else:
            conv = ureg(self._ex_amplitude_unit).to(unit).magnitude
            return self.ex_amplitude_err * conv

    def get_ex_power(self, unit=None):
        if unit is None:
            conv = 1.0
        else:
            conv = ureg(self._ex_power_unit).to(unit).magnitude
        return self.ex_power * conv

    def get_ex_power_err(self, unit=None):
        if unit is None:
            conv = 1.0
        else:
            conv = ureg(self._ex_power_unit).to(unit).magnitude
        return self.ex_power_err * conv

    def reset_masks(self):
        for psd in self.psds.values():
            psd.reset_mask()

    def exclude_freq(self, f_ex, names=None):
        """
        Exclude data points at frequencies f_ex.

        Arguments
        ---------
        f_ex : float or list(floats) or None
            Frequencies to be excluded. If None, the data point at the
            excitation frequency of the excited axis is excluded.
        name : str
            Name of the psd where the data point shall be excluded. If None,
            all psds get f_ex excluded
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            self.psds[name].exclude_freq(f_ex)

    def exclude_freq_outside(self, fmin, fmax, names=None, reset_mask=False):
        """
        Exclude values outside the range fmin, fmax.
        """
        if names and not isinstance(names, list):
            names = [names]

        for name in names if names else self.names:
            if reset_mask:
                self.psds[name].reset_mask()
            self.psds[name].exclude_freq_outside(fmin, fmax)

    def plot_psds(self, names=None, axis=None, **kwargs):
        """
        Plots all power spectral densities or the specified name only.

        Calls the function plot_psd of the PSD object.

        Arguments
        ---------
        name : str or None
        axis : matplotlib.Axis or None
            Axis to add the plot to.
        kwargs : keyword arguments passed to PSD.plot_psd()
            e.g.:
             - plot_all=False
             - plot_masked=False
             - plot_errors=False
             - plot() keyword argumets

        Returns the figure object.
        """
        if names and not isinstance(names, list):
            names = [names]

        if 'title' not in kwargs.keys():
            kwargs['title'] = ('PSDs at {0:1.3f} µm'.format(self.exp_setting.get_height(unit='um')))

        for name in names if names else self.names:
            fig = self.psds[name].plot_psd(axis=axis, **kwargs)
            if axis is None:
                axis = fig.axes[0]

        return fig


def gen_psdm_from_region(region, T_msr, N_avg,
                         T_delay=0.0, psd_traces=None,
                         active_calibration=False, position_traces=None,
                         ex_freq=None,  position_unit='um', exp_setting=None):
    """
    """
    pm = PSDMeasurement(exp_setting=exp_setting, warn=True)

    psd_traces = psd_traces or region.traces
    N_min = int(T_delay * region.samplingrate)
    N_max = int(T_msr * region.samplingrate)
    samples = slice(N_min, N_min + N_max)
    data = region.get_data(psd_traces, samples=samples)

    for name, dat in zip(psd_traces, data.T):
        p = gen_PSD_from_time_series(dat, region.samplingrate, N_avg)
        pm.add_psd(name, p)

    if active_calibration:
        if not ex_freq:
            raise Exception('Excitation frequency must be provided for an'
                            ' active calibration.')

        position_traces = position_traces or ['positionX', 'positionY']
        ex_axis = region._excited(position_traces)
        stage_ex_axis = position_traces[ex_axis]
        stage_signal = region.get_data(stage_ex_axis, samples=samples)[:, 0]
        p = gen_PSD_from_time_series(stage_signal,
                                     region.samplingrate,
                                     N_avg,
                                     calc_errors=True)

        ## !!! important this is already in µm, due to the setup-specific
        ## config file
        ex_amp = float(sqrt(2 * p.psd[p.freq == ex_freq] * p.df))
        ex_amp_err = float(sqrt(2 * p.psd_err[p.freq == ex_freq] * p.df))

        psd_ex_axis = psd_traces[ex_axis]

        freqs = pm.get_freq(psd_ex_axis)
        ex_pow = float(pm.get_psd(psd_ex_axis)[freqs == ex_freq] *
                       pm.psds[psd_ex_axis].df)
        ex_pow_err = float(pm.get_psd_err(psd_ex_axis)[freqs == ex_freq] *
                           pm.psds[psd_ex_axis].df)

        pm.set_ac_params(psd_ex_axis, ex_freq, ex_amp, ex_pow,
                         ex_amplitude_err=ex_amp_err, ex_power_err=ex_pow_err,
                         amplitude_unit=position_unit)

    return pm
#####--------------------------------------------------------------------------
#---- shape ---
#####--------------------------------------------------------------------------

def lorentzian_psd(freq, D, f_c):
    """
    Return the values for frequenzies f of a Lorentzian-shaped **one-sided**
    power spectral density with a diffusion constant D and a corner frequency
    f_c.

    The function applied is: psd = D / (pi**2 * (freq**2 + f_c**2)).

    Arguments
    ---------
    freq : array(float)
        Frequency vector.
    D : float
        Diffusion constant.
    f_c : float
        Corner frequency.

    Returns
    -------
    array

    References
    ----------
    Equ. (9) in:
        Tolić-Nørrelykke et al. (2006)
        Calibration of optical tweezers with positional detection in the back
        focal plane. Review of Scientific Instruments, 77(10), 103101.
        http://doi.org/10.1063/1.2356852
    """
    l = D / (pi**2 * (freq**2 + f_c**2))
    return l


def hydro_psd(freq,  # Hz
              D,  # arb**2/s
              f_c,  # Hz
              radius=0.5e-6,  # m
              height=inf,  # m
              temp=293.15,  # K
              rho=1000,  # kg/m**3
              density_med=None,  # kg/m**3
              viscosity=None,  # Pa*s
              lateral=True,
              verbose=False):
    """
    The hydrodynamically correct power spectral density of a sphere at a given
    height in a viscous medium.

    The function takes Faxén's law and the mass of the sphere into account.

    Arguments
    ---------
    freq : array(float)
        Frequency vector in Hertz
    D : float
        Diffusion coefficient (with no correction, i.e. pure Stokes drag).
    f_c : float
        Corner frequency in Hz.
    radius : float
        Radius of the sphere in meters.
    height : float
        Height, i.e. the bead-center -- surface distance (the real one) in
        meters.
    temp : float
        Absolute temperature in Kelvin
    rho : float
        Mass density of the sphere in kg/m³
    density_med : float
        Mass density of the medium in kg/m³. If None, the density of water at
        the given temperatrue is used.
    viscosity : float
        Viscosity of the medium in Ns/m². If None, the viscosity of water at
        temperature temp is assumed.

    References
    ----------
    Appendix D in:
        Tolić-Nørrelykke et al. (2006)
        Calibration of optical tweezers with positional detection in the back
        focal plane. Review of Scientific Instruments, 77(10), 103101.
        http://doi.org/10.1063/1.2356852

    See Also
    --------
    viscosity_H2O

    """
    drag_stokes = drag(radius, temp, density=density_med,
                       viscosity=viscosity)

    f_c0 = f_c

    drag_l = drag(radius, temp, freq=freq, height=height, density=density_med,
                  viscosity=viscosity, lateral=lateral, verbose=verbose)

    rel_drag = drag_l / drag_stokes

    # f_m0 = drag(radius) / (2 * pi * m)
    # f_m0 = drag(radius) / (2 * pi * rho * (4/3) * pi * radius**3)
    f_m0 = drag_stokes / ((8 / 3) * pi**2 * rho * radius**3)

    P = (D * rel_drag.real /
         (pi**2 * ((f_c0 + freq * rel_drag.imag - freq**2 / f_m0)**2 +
                   (freq * rel_drag.real)**2)))

    return P.real


def low_pass_filter(freq, f3dB, alpha=0):
    """
    Produces the relative PSD of a signal that, to a factor (1-alpha²), is
    filtered by a 1st-order low-pass filter.

    This, in particular, describes the filtering of a quadrant photo-diode:
    Only a fraction of the photons that reach the diode is absorped in the
    depletion region, the other fraction get's absorped outside that region.
    The electron-hole pairs need to diffuse to the depletion region until
    they produce a photo-current. This process causes an effective low-pass
    filter for these photons. The in the power spectrum the function looks
    like this:

     F(f) = alpha² + (1 - alpha²) / (1 + (f / f3dB)**2)

    Arguments
    ---------
    freq : array(float)
        Frequency vector.
    f3dB : float
        Cut-off frequency of the low-pass filter.
    alpha : float
        Filter efficiency. Only a fraction (0 <= alpha <= 1) of the signal that
        is not low pass filtered. Thus, alpha = 1, will produce no low-pass
        filter at all, whereas alpha = 0, will produce a 1st-order low-pass
        filter.

    References
    ----------
    Equ. (35) in:
        Berg-Sørensen, K., & Flyvbjerg, H. (2004)
        Power spectrum analysis for optical tweezers. Review of Scientific
        Instruments, 75(3), 594–612.
        http://doi.org/10.1063/1.1645654
    Equ. (20) in:
        Berg-So̸rensen, K., et al. (2006)
        Power spectrum analysis for optical tweezers. II: Laser wavelength
        dependence of parasitic filtering, and how to achieve high bandwidth.
        Review of Scientific Instruments, 77(6), 063106.
        http://doi.org/10.1063/1.2204589
    """
    return alpha**2 + (1 - alpha**2) / (1 + (freq / f3dB)**2)


def apply_low_pass_filter(fun, f3dB, alpha):
    """
    Decorator function to produces a low-pass filtered modification of the
    input PSD function

    See also
    --------
    low_pass_filter()

    Arguments
    ---------
    fun : function
        Function with first positional argument being the frequency vector
        'freq'.
    f3dB : float
        Cut-off frequency of the low-pass filter.
    alpha : float
        Filter efficiency. Only a fraction (0 <= alpha <= 1) of the signal that
        is not low pass filtered. Thus, alpha = 1, will produce no low-pass
        filter at all, whereas alpha = 0, will produce a 1st-order low-pass
        filter.

    Returns
    -------
    function with arguments (freq, *args, **kwargs).

    References
    ----------
    Equ. (35) in:
        Berg-Sørensen, K., & Flyvbjerg, H. (2004)
        Power spectrum analysis for optical tweezers. Review of Scientific
        Instruments, 75(3), 594–612.
        http://doi.org/10.1063/1.1645654
    Equ. (20) in:
        Berg-So̸rensen, K., et al. (2006)
        Power spectrum analysis for optical tweezers. II: Laser wavelength
        dependence of parasitic filtering, and how to achieve high bandwidth.
        Review of Scientific Instruments, 77(6), 063106.
        http://doi.org/10.1063/1.2204589
    """
    def lp_filtering(freq, *args, **kwargs):
        f = (low_pass_filter(freq, f3dB, alpha=alpha) *
             fun(freq, *args, **kwargs))
        return f

    if hasattr(fun, '__name__'):
        lp_filtering.__name__ = ('low-pass filtered "{0:s}"'
                                 ''.format(fun.__name__))
    if hasattr(fun, '__doc__'):
        lp_filtering.__doc__ = ('Low pass filtered function with cut-off '
                                'frequency f3dB = {0:1.1f} Hz and efficiency '
                                'alpha = {1:1.3f} \n\n'
                                'calls\n{2}{3} \n\n'
                                'Documentation of {2}:\n{4}'
                                ''.format(f3dB, alpha, fun.__name__,
                                          str(signature(fun)), fun.__doc__))
    return lp_filtering


def apply_aliasing(fun, f_sample, N=9):
    """
    Decorator function to produce an aliased version of the input PSD function.

    The aliased version of a psd is calculated by adding the psd values of
    frequencies that lie beyound the sampling frequency to the range below the
    sampling frequency.

     PSD_aliased(f) = sum(PSD(f + n * f_sample)) from n=-N to N.

    Note
    ----
    Note, that N is actually infinity, but a fintite number of summations
    should be sufficient to account for aliasing. N=5 gives about 1.5% error at
    high frequencies, N=10 should give less than 0.5%.

    Arguments
    ---------
    fun : function
        Function to be aliased.
    f_sample : float
        Sampling frequency
    N : int
        Number that defines how many ranges of f_sample should be taken into
        account. The default N=9 give a very good approximation with deviations
        of less than 0.5%
    *args
        passed to the given function.
    **kwargs
        passed to the given function.

    Returns
    -------
    function with arguments (freq, *args, **kwargs)

    References
    ----------
    Equ. (37) in:
        Berg-Sørensen, K., & Flyvbjerg, H. (2004)
        Power spectrum analysis for optical tweezers. Review of Scientific
        Instruments, 75(3), 594–612.
        http://doi.org/10.1063/1.1645654
    """
    def aliasing(freq, *args, **kwargs):
        aliased = sum(fun(freq + i * f_sample, *args, **kwargs)
                      for i in range(-N, N+1))
        return aliased

    if hasattr(fun, '__name__'):
        aliasing.__name__ = 'aliased "{0:s}"'.format(fun.__name__)
    if hasattr(fun, '__doc'):
        aliasing.__doc__ = ('Aliased function with\n'
                            'N = {0} and f_sample = {1:1.1f} Hz. \n\n'
                            'Calls\n'
                            '{2}{3}\n\n'
                            'Documentation of {2}:\n{4}'
                            ''.format(N, f_sample, fun.__name__,
                                      str(signature(fun)), fun.__doc__))
    return aliasing
