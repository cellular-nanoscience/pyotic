# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:29:19 2016

@author: Tobias Jachowski
"""
import numpy as np
import persistent

from .. import helpers as hp


class CalibrationSource(persistent.Persistent):
    """
    Offers attributes to be used by a Calibration.

    Attributes
    ----------
    corrfactor : float
        The correction factor of the size of the bead (i.e. the relative
        stokes drag coefficient far in the solution).
    dsurf : float
        The value of positionZ, where the center of the bead would be on the
        surface.
    beta : numpy.ndarray(n)
        The displacement sensitivity at the position, where the center of the
        bead would be on the surface (nm/mV).
    kappa : numpy.ndarray(n)
        The stiffness at the position, where the center of the bead would be on
        the surface (pN/nm).
    mbeta : numpy.ndarray(n)
        The slope of the displacement sensitivity in units of
        beadcenter-surface distance (i.e. focalshift corrected units of
        positionZ, nm/mV/µm).
    mkappa : numpy.ndarray(n)
        The slope of the stiffness in units of beadcenter-surface distance
        (i.e. focalshift corrected units of positionZ).
    radiusspec : float
        The radius, as it was specified by the manufacturer (µm).
    focalshift : float
        The focalshift of the Setup.
    """
    name = "Standard calibration source"

    def __init__(self, beta=None, kappa=None, radiusspec=0.0, corrfactor=1.0,
                 focalshift=1.0, mbeta=None, mkappa=None, dsurf=0.0, name=None,
                 **kwargs):
        """
        Parameters
        ----------
        beta : tuple or array (n)
            The displacement sensitivity at the position, where the center of
            the bead would be on the surface (nm/mV).
        kappa : tuple or array (n)
            The stiffness at the position, where the center of the bead would
            be on the surface (pN/nm).
        radiusspec : float
            The radius, as it was specified by the manufacturer (µm).
        focalshift : float
            The focalshift of the Setup.
        mbeta : tuple or array (n), optional
            The slope of the displacement sensitivity in units of
            beadcenter-surface distance (i.e. focalshift corrected units of
            positionZ, nm/mV/µm), (numpy.zeros(n), default).
        mkappa : tuple or array (n), optional
            The slope of the stiffness in units of beadcenter-surface distance
            (i.e. focalshift corrected units of positionZ, pN/nm/µm),
            (numpy.zeros(n), default).
        corrfactor : float, optional
            The correction factor of the size of the bead (i.e. the relative
            stokes drag coefficient far in the solution), (1.0, default).
        dsurf : float, optional
            The value of positionZ, where the center of the bead would have
            been on the surface during the calibration measurement (0.0,
            default).
        """
        # Set beta and kappa to generice values. The calibration factors are
        # choosen in such a way that modules relying on a calibration are still
        # (somewhat) functional.
        if beta is None:
            print("Created generic values for beta. Make sure to provide a "
                  "proper calibration, if needed.")
            beta = [1.0, 1.0, 1.0]  # nm/mV
        if kappa is None:
            print("Created generic values for kappa. Make sure to provide a "
                  "proper calibration, if needed.")
            kappa = [1.0, 1.0, 1.0]  # pN/nm
        self.beta = np.array(hp.listify(beta))
        self.kappa = np.array(hp.listify(kappa))
        if mbeta is None:
            mbeta = [0.0 for i in range(self.beta.size)]
        self.mbeta = np.array(hp.listify(mbeta))
        if mkappa is None:
            mkappa = [0.0 for i in range(self.kappa.size)]
        self.mkappa = np.array(hp.listify(mkappa))
        self.dsurf = dsurf
        self.radiusspec = radiusspec
        self.corrfactor = corrfactor
        self.focalshift = focalshift
        self.name = name or ''
