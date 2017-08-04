# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:29:39 2016

@author: Tobias Jachowski
"""
import numpy as np
from scipy import constants
from collections import namedtuple


ForceExtension = namedtuple('ForceExtension', ['extension', 'force'])


def force_extension(bps=3615, T=298.2, samples=1000, min_ext=0.5,
                    max_ext=0.978, pitch=0.338e-9, L_p=43.3e-9, K_0=1246e-12):

    # 0.5-0.978 entspricht ~520nm-1057nm (0.118...49.2pN)
    z = np.linspace(min_ext, max_ext, samples)

    # bps = 3615
    # bps = 1399
    # T=298.2  # K temperature

    k = constants.value('Boltzmann constant')
    # pitch=0.338e-9  # m/bp length of one bp (Saenger 1988)
    # L_p=43.3e-9  # m persistence length
    L_0 = pitch*bps  # m contour length
    # K_0=1246e-12  # N elastic modulus

    F = (k * T / L_p) * (1 / (4 * (1 - z)**2) - 1/4 + z)
    # z = x/L_0 -F/K_0

    x = (z + F / K_0) * L_0

    x = x*1e9  # nm
    F = F*1e12  # pN

    return ForceExtension(extension=x, force=F)


def twistable_wlc(bps=1000, T=298.2, samples=1000, min_F=10e-12, max_F=60e-12,
                  pitch=0.338e-9, L_p=43.3e-9, K_0=1246e-12, S=1500e-12,
                  C=440e-30):
    """
    Twistable worm like chain model Gross et al. 2011
    """
    L_0 = pitch * bps  # m contour length
    k = constants.value('Boltzmann constant')

    # x = np.linspace(min_ext, max_ext, samples)
    F = np.linspace(min_F, max_F, samples)

    def g(F):
        return (S * C - C * F * (x/L_0 - 1 + 1/2
                * (k*T/(F*L_p))**(1/2))**(-1))**(1/2)

    x = L_0 * ((1 - 1/2) * (k * T / (F * L_p))**(1/2)
               + C / (-g(F)**2 + S * C) * F)

    return x
