# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:12:08 2016

@author: Tobias Jachowski
"""
import importlib
import numbers
import numpy as np
import persistent

from . import calibsource as cs
from .. import config as cf
from .. import helpers as hp
from ..graph import GraphMember


class Calibration(GraphMember, persistent.Persistent):
    """
    A calibration gives you easy acces to the stiffness (kappa) and
    displacement sensitivity (beta) at the position, where the microsphere
    would touch the glass surface (intercept) and the dependency of kappa and
    beta of the position in Z (slope).

    A Calibration consists of:
      - a CalibrationSource, containing the distance dependent calibration
        factors
      - the positionZ dependent linear function of stiffness and displacement
        sensitivity
    """
    def __init__(self, calibsource, **kwargs):
        super().__init__(max_parents=0, **kwargs)

        if calibsource is None:
            raise TypeError("Calibration missing the required positional "
                            "argument 'calibsource'.")
        self.calibsource = calibsource

        # If the true distance of the position to the bead changed from the
        # measurement to the calibration (by dh), due to thermal expansion of
        # the stage/objective, the touchdown determined by the touchdown
        # modification corresponds to the dsurf determined by the calibration
        # script and therefore, no offset modification to the linear function
        # has to be applied:
        # h_stage - h_surf - dh != h_stage - h_surf - touchdown
        # If, however, there was no thermal expansion, then the linear function
        # (the intercept) has to be corrected for the touchdown offset:
        # beta0 = beta0 + mbeta*touchdown; kappa0 = kappa0 + mkappa*touchdown
        # To check, whether thermal expansion took place, one would need to
        # save the sum signal of the QPD during calibration and compare/map it
        # to the one saved during the measurement.Percent error of accidently
        # neglected
        # touchdown offset modification for the intercept:
        # offset != touchdown
        # self.offset = offset  # um, in case dsurf was not determined
        # correctly: offset of self.calibsource.dsurf, is set according to
        # touchdown of class Touchdown
        # self.correct_offset = True  # formerly thermal_expansion !=
        # correct_offset
        # correct the radius or assume the specified radius to be the real one?
        self.correct_radius = True

    def beta(self, positionZ=0.0):
        """
        Return the positionZ dependent displacement sensitivity.
        """
        return self._beta[0] + positionZ * self._beta[1]

    def kappa(self, positionZ=0.0):
        """
        Return the positionZ dependent stiffness.
        """
        return self._kappa[0] + positionZ * self._kappa[1]

    @property
    def _beta(self):
        """
        Give the slope (m) and intercept (y_0) of the positionZ dependent
        displacement sensitivity.
        See also intercept and slope for explanation.
        In (nm/mV)
        """
        return np.array([self.intercept('beta'), self.slope('beta')])

    @property
    def _kappa(self):
        """
        Give the slope (m) and intercept (y_0) of the positionZ dependent
        stiffness.
        See also intercept and slope for explanation.
        In (pN/nm)
        """
        return np.array([self.intercept('kappa'), self.slope('kappa')])

    def intercept(self, para):
        """
        Get the displacement sensitivity (para='beta') or stiffness
        (para='kappa') at the position, where the bead surface would touch the
        glass surface (the intercept). The position and in turn the value
        depends on the radius of the bead and therefore, whether the corrected
        radius or the specified radius is used for the calibration (choosen by
        attribute self.correct_radius).
        """
        # Give the intercept `y_0` for the function of 'surface distance' and
        # 'bead radius offset' corrected 'positionZ' (um)
        # (bead-surface glass-surface) for displacement sensitivity
        # (para='beta') and stiffness (para='kappa').
        # f(h_stage - dsurf - radius/focalshift) = y_0 + m*d
        if para == 'beta':
            return (self.calibsource.beta + self.calibsource.mbeta
                    * self.radius)
        elif para == 'kappa':
            return (self.calibsource.kappa + self.calibsource.mkappa
                    * self.radius)

    def slope(self, para):
        """
        Get the positionZ dependent slope of the displacement sensitivity
        (para='beta') or stiffness (para='kappa').
        """
        # Gives the slope `m` for the function of 'surface distance' and
        # 'bead radius offset' corrected 'positionZ' (um)
        # (bead-surface to glass-surface) for displacement sensitivity
        # (para='beta') and stiffness (para='kappa').
        # f(h_stage - dsurf - radius/focalshift) = y_0 + m*d
        # positionZ ~ - distance of calibsource -> negative slope
        if para == 'beta':
            return - (self.calibsource.mbeta * self.calibsource.focalshift)
        elif para == 'kappa':
            return - (self.calibsource.mkappa * self.calibsource.focalshift)

    def displacement(self, psd, positionZ=0.0, calib_traces=None):
        """
        Calculate the displacement.

        The displacement of psd is calculated with the displacement
        sensitivity stored in self.beta(). Height dependency of the
        displacement sensitivity is accounted for by supplying a positionZ
        signal.

        Parameters
        ----------
        psd : 2D array
            A 2D array containing the psd data.
        positionZ : int or array, optional
            Integer or array containing the height where each value of psd was
            measured. If an array is supplied, it needs to have the same length
            as psd.
        calib_traces : int, list, or slice, optional
            If self.beta() has another length or another order than the second
            dimension of psd, than you can specify which psd traces correspond
            to which self.beta() values. calib_traces needs to have as many
            values, as psd.shape[1].

        Returns
        -------
        2D array
            An array with the same dimensions as psd with displacement values.
        """
        if calib_traces is None:
            calib_traces = slice(0, psd.shape[1])
        if isinstance(positionZ, numbers.Number):
            positionZ = np.array([[positionZ]])
        if positionZ.ndim < 2:
            positionZ = positionZ[:, np.newaxis]
        beta = self.beta(positionZ=positionZ)[:, calib_traces]
        displacement = beta * psd
        return displacement

    def force(self, displacement, positionZ=0.0, calib_traces=None):
        """
        Calculate the force.

        The force of displacement is calculated with the stiffness stored in
        self.kappa(). Height dependency of the stiffness is accounted for by
        supplying a positionZ signal.

        Parameters
        ----------
        displacement : 2D array
            A 2D array containing the displacement (see self.displacement()).
        positionZ : int or array
            Integer or array containing the height where each value of
            displacement was measured. If an array is supplied, it needs to
            have the same length as psd.
        calib_traces : int, list, or slice
            If self.kappa() has another length or another order than the second
            dimension of displacement, than you can specify which displacement
            traces correspond to which self.kappa() values. calib_traces needs
            to have as many values, as psd.shape[1].

        Returns
        -------
        2d array
            An array with the same dimensions as displacement with force
            values.
        """
        if calib_traces is None:
            calib_traces = slice(0, displacement.shape[1])
        if isinstance(positionZ, numbers.Number):
            positionZ = np.array([[positionZ]])
        if positionZ.ndim < 2:
            positionZ = positionZ[:, np.newaxis]
        kappa = self.kappa(positionZ=positionZ)[:, calib_traces]
        force = kappa * displacement
        return force

    @property
    def dsurf(self):
        """
        Return the distance, where the bead center would have been on the glass
        surface during the calibration measurement.
        """
        return self.calibsource.dsurf

    @property
    def psurf(self):
        """
        Return the positionZ, were the bead center would have been on the glass
        surface during the calibration measurement.
        """
        return - self.dsurf

    @property
    def touchdown(self):
        """
        Return the positionZ, where the bead surface would have touched the
        glass surface during the calibration measurement.
        Uses either the corrected or specified radius of the bead determined by
        attribute `self.correct_radius`.
        """
        # Keep in mind, dsurf from the calibration script relates to distance,
        # which is ~ to - positionZ
        # -> !-! (dsurf !+! self.radius/self.focalshift)
        # = psurf - self.radius/self.focalshift
        return self.psurf - self.radius / self.focalshift

    @property
    def correct_radius(self):
        return self._correct_radius

    @correct_radius.setter
    def correct_radius(self, value):
        if not hasattr(self, '_correct_radius') \
                or value != self._correct_radius:
            self._correct_radius = value
            self.set_changed()

    @property
    def radius(self):
        if self.correct_radius:
            radius = self.calibsource.radiusspec * self.calibsource.corrfactor
        else:
            radius = self.calibsource.radiusspec
        return radius

    @property
    def focalshift(self):
        return self.calibsource.focalshift


def create_calibration(source_type=None, source_module=None, source_class=None,
                       cfgfile=None, verbose=False, name=None, **kwargs):
    """
    Parameters
    ----------
    source_type : str, optional
        Choose between different types of calibrationsources to be created or
        loaded. If None is selected and module or source_class are not given, a
        generic calibrationsource is created. The calibration factors are
        choosen in such a way that modules relying on a calibration are still
        (somewhat) functional.
        If a calibrationsource type is choosen, which needs to read in a file,
        additionally the parameters 'directory' and 'filename' have to be
        provided.
    source_module : str, optional
        A str of the name of a module a source_class is located in.
    source_class : str, optional
        A str of a class tha should be used as a CalibrationSource. The class
        has to be located in the module described.
    **kwargs : optional
        Any parameters, a CalibrationSource needs to be initialized.
    """
    # Set default configfile
    cfgfile = cfgfile or 'calibration.cfg'

    # Read configfile
    cfg = cf.read_cfg_file(cfgfile)

    # If directory and/or filename are given, join them to full path filename
    filename = kwargs.pop('filename', None)
    directory = kwargs.pop('directory', None)
    if filename is not None:
        _f, _d, filename = hp.file_and_dir(filename, directory)
        kwargs['filename'] = filename

    # Create a cs_class. If no source_module and source_class are given
    # by user, default to configfile. If configfile doesn't have entries
    # for source_module and source_class, default to std_mod and std_cls
    if source_module is None or source_class is None:
        # Get calibsource class
        cs_class = source_class or cf.get_cfg_class(
            cfg,
            sec=source_type,
            mod_opt='source_module',
            cls_opt='source_class',
            std_mod='.calibration',
            std_cls='CalibrationSource',
            verbose=verbose)
    else:
        mod = importlib.import_module(source_module, __package__)
        cs_class = getattr(mod, source_class)

    # Get parameters for cs_class
    _kwargs = cf.get_cfg_sec_dict(cfg, sec=source_type, verbose=verbose)

    # Remove options for cs_class creation (module and class)
    _kwargs.pop('source_module', None)
    _kwargs.pop('source_class', None)

    # Convert all comma separated strings into list of floats
    for key, value in _kwargs.items():
        if isinstance(value, str):
            _kwargs[key] = [float(l.strip()) for l in value.split(',')]

    # update options for cs_class, user given kwargs take precedence over
    # options retrieved from configfile
    _kwargs.update(kwargs)

    # If it is a standard calibration, put the source_type into the name.
    cs_name = None
    if cs_class is cs.CalibrationSource and source_type is not None:
        cs_name = ''.join(('Standard calibration source: ', source_type))

    # Create calibrationsource
    calibsource = cs_class(name=cs_name, **_kwargs)

    # remove parameters used by calibsource
    for par in ['beta', 'kappa', 'radiusspec', 'focalshift',
                'mbeta', 'mkappa', 'corrfactor', 'dsurf']:
        _kwargs.pop(par, None)

    # create calibration
    return Calibration(calibsource, name=name, **_kwargs)
