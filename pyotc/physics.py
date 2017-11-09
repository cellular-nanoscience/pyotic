# -*- coding: utf-8 -*-
# """
# - Author: steve simmert
# - E-mail: steve.simmert@uni-tuebingen.de
# - Copyright: 2015
# """
"""
Provides physical functions and values.
"""
from . import ureg
from .utilities import str2u

from scipy import complex_
from scipy import exp
from scipy import float_
from scipy import inf
from scipy import isinf
from scipy import pi
from scipy import sqrt

import warnings


MATERIALS = {'PS': 1050,  # kg/m³
             'Polystyrene': 1050,
             'polystyrene': 1050,
             'ps': 1050,
             'Silica': 2000,
             'silica': 2000,
             'SiO': 2000,
             'sio': 2000,
             'Titania': 4230,
             'titania': 4230,
             'TiO': 4230,
             'tio': 4230}


class Material(object):
    """
    Material of the particle.
    """
    def __init__(self, name, density=None, density_unit='kg/m**3'):
        """
        Specify the name of the material of the particle.

        If the name is unknown you need to specify the density, as well.
        """
        self.name = name
        if name in MATERIALS:
            self._density = MATERIALS[name]
            self._density_unit = 'kg/m**3'
        else:
            if density:
                self._density = density
                self._density_unit = str2u(density_unit)
            else:
                raise Exception('Unknown Material "{}" and density is not'
                                'given.'.format(name))

    def get_density(self, unit=None):
        if unit is None or unit == self._density_unit:
            conv = 1.0
        else:
            conv = ureg(self._density_unit).to(unit).magnitude
        return self._density * conv


def viscosity_H2O(temp):
    """
    Calculates the viscosity (Pa*s) of water at temperature T in Kelvin.

    Warning
    -------
    Function is a fit to data of the viscosity of water at different
    temperatures [0, 100].

    Do not use it outside this range!

    References
    ----------
    CRC Handbook of Chemistry and Physics
        88th Edition
        CRC Press; 88 edition (October 1, 2007)
    """
    a = 3.4175e-4
    b = -0.01276
    c = 3.5857e-5

    if not(273.15 <= temp <= 373.15):
        warnings.warn('This function is only valid within a temperature'
                      'range between [273.15, 373.15] Kelvin!\n'
                      'T = {0:1.3f}'.format(temp))

    return a / (1 + b * temp + c * temp**2)


def dviscosity_H2O(temp):
    """
    Calculates the slope of the viscosity of water at a given temperature.

    See also
    --------
    viscosity_H2O
    """
    # the values to calculate vsicosity in Pa*s from T in K
    a = 3.4175e-4
    b = -0.01276
    c = 3.5857e-5

    if not(273.15 <= temp <= 373.15):
        warnings.warn('This function is only valid within a temperature'
                      'range between [273.15, 373.15] Kelvin!\n'
                      'T = {0:1.3f}'.format(temp))

    out = -((a * (b + 2 * c * temp)) / (1 + b * temp + c * temp**2)**2)
    return out


def density_H2O(temp):
    """
    Returns the density of water in kg/m³ at the specified temperature in
    Kelvin.

    Warning
    -------
    This function relies on a 3rd-order polynomial-fit to experimental
    data within the temperature region [0, 100] degrees Celsius. Within that
    range the fit reproduces the data within < 0.03 percent!

    Do not use it outside this range!

    References
    ----------
    CRC Handbook of Chemistry and Physics
        88th Edition
        CRC Press; 88 edition (October 1, 2007)
    """
    T = temp - 273.15  # convert to celcius scale

    if not(0 <= T <= 100):
        warnings.warn('This function is only valid within a tempearture'
                      'range between [273.15, 373,15] Kelvin!\n'
                      'T = {0:1.3f}'.format(T))

    a = 1.5690e-05
    b = -0.00591743
    c = 0.02020695
    d = 999.949648

    rho = a * T**3 + b * T**2 + c * T + d

    return rho


def kinematic_viscosity_H2O(temp):
    """
    Return the kinematic viscosity of water (m²/s) at the specified
    temperature in Kelvin.

    Warning
    -------
    Note that this function relies on fits to experimental
    data within the temperature region [0, 100] degrees Celsius. Within that
    range the fit reproduces the data within < 0.03 percent!

    Do not use it outside this range!

    References
    ----------
    CRC Handbook of Chemistry and Physics
        88th Edition
        CRC Press; 88 edition (October 1, 2007)
    """
    if not(0 <= temp-273.15 <= 100):
        warnings.warn('This function is only valid within a tempearture'
                      'range between [273.15, 373,15] Kelvin!')
    eta = viscosity_H2O(temp)
    rho = density_H2O(temp)

    return eta / rho


def drag(radius, temp, freq=0.0, height=inf, density=None, viscosity=None,
         lateral=True, verbose=False):
    """
    Calculate the frequency-dependent viscous drag of a sphere close to an
    infinit plane surface.

    The function uses an approximation to Faxen's solution for the drag near a
    plane surface at frequency f for a sphere of radius R. The drag depends on
    the height (surface to bead-center distance) of the sphere above the
    surface. Further, it depends on the kinematic viscosity of the medium, thus
    you can provide density and viscosity of the medium. If not provided, the
    kinematic viscosity of water at the given temperature is calculated.

    Arguments
    ---------
    radius : float
        Radius of the sphere in meter.
    temp : float
        Absolute temperature in kelvin.
    height : float
        Surface to bead-center distance in meter. If height = inf, the function
        returns stokes drag.
    density : float or None
        Density of the medium in kg/m³. If None, the density of water at the
        given temperature is used.
    viscosity : float
        Viscosity of the medium in Pa*s.  If None, the viscosity of water at
        the given temperature is used.
    later : bool
        !Not implemented yet!
        If true, the lateral drag is calculated. If false, the axial one
        is calculated.

    [Ref]: Tolić-Nørrelykke et al. Rev. Sci. Instrum. 77, 103101 (2006)

    Nota bene:
    In contrast to [Ref.], here f_nu is directly calculated with the given
    viscosity.
    """
    R = float_(radius)
    T = float_(temp)
    f = float_(freq)
    l = float_(height) if height >= radius else radius

    if verbose and height < radius:
        warnings.warn('Height is less than specified radius: '
                      '{0:1.3e} < {1:1.3e}\n'
                      'Height is set to radius as fallback.'
                      ''.format(height, radius))

    # viscous drag far away from a surface at constant speed
    eta = viscosity if viscosity else viscosity_H2O(T)
    d_0 = 6 * pi * eta * R

    if density is None or viscosity is None:
        nu = kinematic_viscosity_H2O(T)
        if verbose:
            print('Kinematic viscosity of water is used: '
                  '{0:1.4e} m²/s'.format(nu))
    else:
        nu = viscosity / density
        if verbose:
            print('Given viscosity and density is used '
                  '{0:1.4e} Pa*s / {1:1.4e} kg/m³ = {2:1.4e} m²/s'
                  ''.format(viscosity, density, nu))

    # f_nu - characteristic frequency
    f_nu = nu / (pi * R * R)

    # relative frequency
    f_ = f / f_nu

    # viscous drag far away from a surface at sinusodial movement
    # at frequency f
    # [Ref.] Equ. D4
    d_stokes = d_0 * (1 + complex_(1-1j) * sqrt(f_) - complex_(1j) * (2/9) * f_)

    if isinf(l):
        return d_stokes

    if lateral:
        # viscous drag at distance l from the surface at sinusodial movement f
        # [Ref.] Equ. D6
        delta = R / sqrt(f_)

        expon = 1 - exp((-(2 * l - R) / delta) * complex_(1-1j))

        # workaround 1 / sqrt(complex_(0)) = inf + nan j
        try:
            expon[f_ == 0.0] = complex_(0)
        except:
            if f_ == 0.0:
                # assume f_ is a scalar
                expon = complex_(0)

        denom = (1 - (9/16) * (R/l) * (1 - sqrt(f_)/3 * complex_(1-1j) +
                                       f_ * complex_(2j/9) - (4/3) * expon))

        d = d_stokes / denom
    else:
        # viscous drag at distance l from the surface at sinusodial movement f
        # [Ref.] Equ. D6
        delta = R / sqrt(f_)

        expon = 1 - exp((-(2 * l - R) / delta) * complex_(1-1j))

        # workaround 1 / sqrt(complex_(0)) = inf + nan j
        try:
            expon[f_ == 0.0] = complex_(0)
        except:
            if f_ == 0.0:
                # assume f_ is a scalar
                expon = complex_(0)

        denom = (1 - (9/16) * (R/l) * (1 - sqrt(f_)/3 * complex_(1-1j) +
                                       f_ * complex_(2j/9) - (4/3) * expon))

        d = d_stokes / denom

    return d


def faxen_factor(height, radius, lateral=True):
    """
    Return the relative change of the stokes drag according to Faxen's law.

    References
    ----------
    Schäffer, E., et al.
        Surface forces and drag coefficients of microspheres near a plane
        surface measured with optical tweezers.
        Langmuir 23, 3654–3665 (2007).
    """
    h = height
    r = radius

    if lateral:
        # Formular (5) in Schäffer et al.
        factor = 1 / (1 - 9 / 16 * r / h
                      + 1 / 8 * (r / h)**3
                      - 45 / 256 * (r / h)**4
                      - 1 / 16 * (r / h)**5
                      )
    else:
        # Formular (6) in  Schäffer et al.
        factor = 1 / (1 - 9 / 8 * r / h
                      + 1 / 2 * (r / h)**3
                      - 57 / 100 * (r / h)**4
                      + 1 / 5 * (r / h)**5
                      + 7 / 200 * (r / h)**11
                      - 1 / 25 * (r / h)**12
                      )
    return factor


def oseen_factor(height, radius, distance):
    """
    Factor that corrects stokes drag for two walls with given distance.

    Arguments
    ---------
    height : float
        height above the lover surface
    radius : float
        radius of the particle
    distance : float
        distance between the two walls

    References
    ----------
    Dufresne, E. R. et al.
        Brownian dynamics of a sphere between parallel walls.
        Europhys. Lett. 53, 264–270 (2001).
        Equation 5
    """
    o = (faxen_factor(height, radius) +
         faxen_factor(distance - height, radius) - 1)
    return o
