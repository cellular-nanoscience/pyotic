# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 13:36:39 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np

from collections.abc import Iterable

from . import signal as sn
from .. import helpers as hp
from .evaluate import Evaluator
from .signalfeature import CycleSectioner


class Tether(Evaluator):
    """
    Calculates tether related parameters, like sections of stressing or
    releasing the bead, bead-center surface distance, extension of the tether,
    or force acting on the tether.

    To distinguish between the excited axes ('x', 'y'), directions
    ('left', 'right'), and different cycles ('stress', 'release'), a
    CycleSectioner is used.
    """
    def __init__(self, region=None, calibration=None, resolution=None,
                 filter_time=-1.0, resolution_sf=None, filter_time_sf=None,
                 **kwargs):
        """
        This __init__() constructor extends the superclass (Evaluator)
        constructor by setting the traces_sf to default to 'positionXY'.

        For extensive documentation of the Parameters and additional
        attributes, see `pyoti.evaluate.evaluate.Evaluator`.

        Parameters
        ----------
        region : pyoti.region.region.Region
        calibration : pyoti.calibration.calibration.Calibration
        resolution : float
        filter_time : float
        resolution_sf : float
        filter_time_sf : float
        **kwargs
            Extra arguments, e.g. traces_sf, in case the positionXY signals are
            called differently.

        Attributes
        ----------
        displacementXYZ : 2D numpy.ndarray
        forceXYZ : 2D numpy.ndarray
        force
        distanceXYZ
        distance
        extension
        angle
        rightstress
        leftstress
        leftrelease
        rightrelease
        stress
        release
        rightstresses
        leftstresses
        leftreleases
        rightreleases
        stresses
        releases
        """
        sf_class = CycleSectioner
        traces_sf = kwargs.pop('traces_sf', 'positionXY')

        # resolution for (almost) all properties, e.g. displacement, force, ...
        # resolution 1.000 Hz, samplingrate 40.000 Hz -> 40 points
        super().__init__(region=region, calibration=calibration,
                         resolution=resolution, filter_time=filter_time,
                         sf_class=sf_class, traces_sf=traces_sf,
                         resolution_sf=resolution_sf,
                         filter_time_sf=filter_time_sf)


    def _sections(self, axis=None, direction=None, cycle=None,
                  concat_direction=False, concat_cycle=False, info=False,
                  **kwargs):
        """
        Calculate start/stop indices (segments) of the requested axis,
        direction, and cycle (sections).

        Parameters
        ----------
        axis : str or list of str, optional
            'x', 'y', or ['x', 'y'] (default)
        direction : str or list of str, optional
            'left', 'right', or ['left', 'right'] (default)
        cycle : str or list of str, optional
            'stress', 'release', or ['stress', 'release'] (default)
        concat_direction : bool, optional
            Concatenate individual left/right sections into one section. Only
            evaluated, if `concat_cycle` is True.
        concat_cycle : bool, optional
            Concatenate individual stress/release sections into one new section
            with the start beeing the start of stress and the stop beeing the
            stop of release.
        info : bool, optional
            Additionally to the sections return a 2D numpy.ndarray, containing
            the axis, direction, and the cycle for every individual section.
            Axis can be either 'x' or 'y'. Direction can be 'left', 'right' or
            'leftright'. Cycle can be 'stress', 'release', or 'stressrelease'.

        Returns
        -------
        2D np.ndarray or tuple of two 2D np.ndarray
            If info is False start/stop values in rows.
            If info is True start/stop values in rows and an 2D np.ndarray with
            infos (see parameter `info`).
        """
        if axis is None:
            axis = ['x', 'y']
        if direction is None:
            direction = ['left', 'right']
        else:
            # A direction was chosen, prevent concat
            concat_direction = False
        if isinstance(direction, str):
            direction = [direction]
        if cycle is None:
            cycle = ['stress', 'release']
        else:
            # A cycle was chosen, prevent concat
            concat_cycle = False
        if isinstance(cycle, str):
            cycle = [cycle]

        sections = np.empty((0, 2), dtype=int)
        infos = np.empty((0, 3), dtype=str)

        # Get all requested segments
        # for ax, axi in zip(axis, range(len(axis))):
        for ax in axis:
            # separate release/stress and separate left/right -> x/y,
            # right/left, stress/release
            if not concat_cycle:
                for sec, dc in zip([d + c for d in direction for c in cycle],
                                   [[d, c] for d in direction for c in cycle]):
                    for segment in self._sf.sections[ax]:
                        sections = np.r_[sections, segment[sec]]
                        if len(segment[sec]) > 0:
                            seg_info = np.array([[ax, dc[0], dc[1]]]
                                                * len(segment[sec]))
                            infos = np.r_[infos, seg_info]
            # concat release/stress and separate left/right -> x/y, right/left
            if concat_cycle and not concat_direction:
                # for d, di in zip(direction, range(len(direction))):
                for d in direction:
                    for segment in self._sf.sections[ax]:
                        sections = np.r_[sections, segment[d]]
                        if len(segment[d]) > 0:
                            seg_info = np.array([[ax, d, 'stressrelease']]
                                                * len(segment[d]))
                            infos = np.r_[infos, seg_info]
            # concat release/stress and concat left/right -> x/y
            if concat_cycle and concat_direction:
                sections = np.r_[sections, self._sf.excited[ax]]
                if len(self._sf.excited[ax]) > 0:
                    seg_info = np.array([[ax, 'leftright', 'stressrelease']]
                                        * len(self._sf.excited[ax]))
                    infos = np.r_[infos, seg_info]

        if info:
            return sections, infos
        else:
            return sections

    def _extrema(self, axis=None, extremum=None):
        if axis is None:
            axis = ['x', 'y']
        if extremum is None:
            extremum = ['minima', 'maxima']
        elif isinstance(extremum, str):
            extremum = [extremum]

        extrema = np.empty(0, dtype=int)

        # Get all extrema and sort
        for ax in axis:
            for section in self._sf.sections[ax]:
                for ext_type in extremum:
                    idx = section[ext_type]
                    idx = self.undecimate_and_limit(idx)
                    extrema = np.r_[extrema, idx]
        extrema.sort()
        return extrema


    def stress_release_pairs(self, i=None, axis=None, direction=None,
                             slices=True, decimate=None, reduce_list=False,
                             **kwargs):
        """
        Calculate start/stop indices (as segments or slices) of stress/release
        cycle pairs of the requested axis and direction.

        Parameters
        ----------
        i : int, optional
            The index of the stress release pair to be selected and returned
            from all stress release pairs calculated. Depending on `axis` and
            `direction` the total availabe number of stress release pairs may
            vary.
        axis : str or list of str, optional
            'x', 'y', or ['x', 'y'] (default)
        direction : str or list of str, optional
            'left', 'right', or ['left', 'right'] (default)
        slices : bool, optional
            Set to False to return segments instead of slices. Default is True.
        decimate : int, optional
            Used to set the step attribute of the returned slices. Only
            evaluated, if `slices` is True.
        **kwargs : dict, optional
            Only used for compatibility purposes, to be able to call method
            with parameters not defined in method definition.

        Returns
        -------
        dict
            stress : 1D np.ndarray of slices or segments
                idx: 1D np.ndarray of slices or segments
                info: 1D np.ndarray of 1D np.ndarrays of type str
                    The str arrays have the form of [axis, direction, cycle].
            release : 1D np.ndarray of slices or segments
                info: 1D np.ndarray of 1D np.ndarrays of type str
                    The str arrays have the form of [axis, direction, cycle].
        """
        axis = ['x', 'y'] if axis is None else axis
        direction = ['left', 'right'] if direction is None else direction
        direction = [direction] if isinstance(direction, str) else direction

        # Get unpaired stress/release indices and initialize paired
        # dictionaries
        pairs = {}
        _pairs = {}
        for cycle in ['stress', 'release']:
            pairs[cycle] = { 'idx': [], 'info': [] }
            _pairs[cycle] = {}
            _pairs[cycle]['idx'], _pairs[cycle]['info'] \
                = self.sections(axis=axis, direction=direction, cycle=cycle,
                                slices=False, info=True)

        # Group all stress/release cycle pairs according to the extrema
        # A stress release pair corresponds to one extremum, only if the stop
        # of the stress and the start of the release segment equal the
        # extremum. However, either one of the stress or release segment can be
        # missing.
        extrema = self._extrema(axis=axis)
        for ext in extrema:
            # Get pair of stress/release whose stop/start is equal to extremum
            i_str = _pairs['stress']['idx'][:, 1] == ext
            i_rls = _pairs['release']['idx'][:, 0] == ext
            stress = _pairs['stress']['idx'][i_str]
            release = _pairs['release']['idx'][i_rls]
            stress_info = _pairs['stress']['info'][i_str]
            release_info = _pairs['release']['info'][i_rls]

            if stress.size == 0 and release.size > 0:
                # Create stress slice, if only release is valid
                stress = np.array([[ext, ext]])
            if release.size == 0 and stress.size > 0:
                # Create release slice, if only stress is valid
                release = np.array([[ext, ext]])
            if stress.size > 0 and release.size > 0:
                pairs['stress']['idx'].append(stress)
                pairs['release']['idx'].append(release)
                pairs['stress']['info'].append(stress_info)
                pairs['release']['info'].append(release_info)

        # Get the maximum number of stress release pairs, prevent index
        # overflow, and allow for negative indices
        if i is None:
            idx = slice(None)
        else:
            if isinstance(i, Iterable):
                idx = i
            else:
                stop = len(pairs['stress']['idx'])
                if i < 0: i = stop + i
                i = min(max(0, i), stop - 1)
                idx = [i]

        # Select idx and convert lists into arrays and segments into slices
        for cycle in ['stress', 'release']:
            for sort in ['idx', 'info']:
                pairs[cycle][sort] = np.concatenate(pairs[cycle][sort])[idx]
            if slices:
                pairs[cycle]['idx'] = \
                    sn.idx_segments_to_slices(pairs[cycle]['idx'],
                                              decimate=decimate)

        # If only one pair, remove list and select only pair
        if reduce_list and len(pairs['stress']['idx']) == 1:
            for cycle in ['stress', 'release']:
                for sort in ['idx', 'info']:
                    pairs[cycle][sort] = pairs[cycle][sort][0]

        return pairs


    def baseline_idx(self, axis=None, strict=False, extrapolate=False):
        """
        Calculate the indices of the baseline, according to the excited axes,
        i.e. the indices, where the excited axis is zero. This method uses
        the detected sections of `self.region`, instead of searching for
        the value 0.0 in the excited axes, as compared to the function
        `pyoti.evaluate.signal.basline_idx()`.

        Parameters
        ----------
        axis : str or list of str, optional
            'x', 'y', or ['x', 'y'] (default)
        strict : bool
            Take only baseline indices of releases's stops or stresses's
            starts, if a stress directly follows a release section.
        extrapolate : bool
            Take baseline indices, even if they are indexed by only one
            stress or release segment, i.e. indices which mark the start or
            the end of a wave.
        """
        # Get stress/release pairs
        if strict:
            # Get only stress/release pairs, i.e. only segment pairs, whose
            # release's stop and stress's start index equal an extremum.
            pairs = self.stress_release_pairs(axis=axis, slices=False)
            stresses = pairs['stress']['idx']
            releases = pairs['release']['idx']
        else:
            # Get all stress/release segments, even those without a
            # corresponding extremum.
            stresses = self.sections(axis=axis, cycle='stress', slices=False)
            releases = self.sections(axis=axis, cycle='release', slices=False)

        # A baseline point lies exactly beetween a release and a following
        # stress section. Therefore, the baseline point idx equals a stress's
        # start and a release's stop.
        stress_base_idx = stresses[:, 0]
        release_base_idx = releases[:, 1]

        # Remove stresses, whose start, and releases, whose stop is equal to
        # an extremum, which means, the excitation as above or below the
        # baseline.
        extrema = self._extrema(axis=axis)
        stress_base_idx = np.setdiff1d(stress_base_idx, extrema,
                                       assume_unique=True)
        release_base_idx = np.setdiff1d(release_base_idx, extrema,
                                        assume_unique=True)

        base_idx = np.r_[stress_base_idx, release_base_idx]
        if strict or not extrapolate:
            # Take only baseline point indices, that lie exactly between a
            # release and a following stress section, i.e. have two equal
            # entries, one from a release's start and one from a stress's stop.
            # I.e. sort out all indices that did not come from a stress/release
            # pair but instead only a single stress or release segment.
            base_idx.sort()
            base_idx = base_idx[np.r_[base_idx[:-1] == base_idx[1:], False]]
        return np.unique(base_idx)

    def _rfigure(self, legend=True, fig=None, ax=None):
        """
        Plot the overview of detected stress release cycles.
        """
        if fig is None and ax is None:
            fig, ax = plt.subplots()
            suptitle = True
        elif fig is None:
            fig = ax.get_figure()
            suptitle = False
        elif ax is None:
            ax = fig.gca()
            suptitle = False

        ax.grid(True)

        line_rstr = None
        line_rrls = None
        line_lstr = None
        line_lrls = None
        line_minima = None
        line_maxima = None
        t = self.timevector
        for axis, trace in zip('xy', ['positionX', 'positionY']):
            s = self.get_data(traces=trace) * 1e6  # m -> µm
            r_str_rls = self.stress_release_pairs(axis=axis, direction='right')
            l_str_rls = self.stress_release_pairs(axis=axis, direction='left')
            rstr = r_str_rls['stress']['idx']
            lstr = l_str_rls['stress']['idx']
            rrls = r_str_rls['release']['idx']
            lrls = l_str_rls['release']['idx']

            ax.plot(t, s, lw=0.1, ms=2, color='k', alpha=1.0)

            # line_rstr = None
            # line_rrls = None
            # line_lstr = None
            # line_lrls = None
            for rstr, rrls in zip(rstr, rrls):
                line_rstr, = ax.plot(t[rstr], s[rstr], lw=0.4, ms=2, color='m')
                line_rrls, = ax.plot(t[rrls], s[rrls], lw=0.4, ms=2, color='c')
            for lstr, lrls in zip(lstr, lrls):
                line_lstr, = ax.plot(t[lstr], s[lstr], lw=0.4, ms=2, color='g')
                line_lrls, = ax.plot(t[lrls], s[lrls], lw=0.4, ms=2, color='y')

            # line_minima = None
            # line_maxima = None
            for segment in self._sf.sections[axis]:
                minima = self.undecimate_and_limit(segment['minima'])
                maxima = self.undecimate_and_limit(segment['maxima'])
                line_minima, = ax.plot(t[minima], s[minima], '.', ms=5,
                                       color='b')
                line_maxima, = ax.plot(t[maxima], s[maxima], '.', ms=5,
                                       color='r')

        line_excited_x = None
        for x_c in (self.undecimate_and_limit(self._sf.excited['x'])
                    / self.resolution):
            line_excited_x = ax.hlines(0.0, x_c[0], x_c[1], alpha=1,
                                       colors='b', linestyle='solid', lw=1)
            # ax.plot(x_c[0], 0.5, '.k', alpha=1, ms=3)
            # ax.plot(x_c[1], 0.5, '.k', alpha=1, ms=3)
            ax.vlines(x_c[0], -0.01, 0.01, alpha=1, colors='b',
                      linestyle='solid', lw=1)
            ax.vlines(x_c[1], -0.01, 0.01, alpha=1, colors='b',
                      linestyle='solid', lw=1)

        line_excited_y = None
        for y_c in (self.undecimate_and_limit(self._sf.excited['y'])
                    / self.resolution):
            line_excited_y = ax.hlines(0.0, y_c[0], y_c[1], alpha=1,
                                       colors='r', linestyle='solid', lw=1)
            # ax.plot(y_c[0], -0.5, '.k', alpha=1, ms=3)
            # ax.plot(y_c[1], -0.5, '.k', alpha=1, ms=3)
            ax.vlines(y_c[0], -0.01, 0.01, alpha=1, colors='r',
                      linestyle='solid', lw=1)
            ax.vlines(y_c[1], -0.01, 0.01, alpha=1, colors='r',
                      linestyle='solid', lw=1)

        ax.set_xlim((t[0], t[-1]))

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Signal positionX and Y (µm)")
        if suptitle:
            fig.suptitle("Automatically detected excited axis, minima, "
                         "maxima, and sections.")

        if legend:
            if line_minima is not None:
                line_minima.set_label('minima')
            if line_maxima is not None:
                line_maxima.set_label('maxima')
            if line_rstr is not None:
                line_rstr.set_label('rightstress')
            if line_rrls is not None:
                line_rrls.set_label('rightrelease')
            if line_lstr is not None:
                line_lstr.set_label('leftstress')
            if line_lrls is not None:
                line_lrls.set_label('leftrelease')
            if line_excited_x is not None:
                line_excited_x.set_label('excited x')
            if line_excited_y is not None:
                line_excited_y.set_label('excited y')

            ax.legend(loc='upper right')

        return fig


    def force_extension_pairs(self, i=None, axis=None, direction=None,
                              decimate=None, twoD=False, posmin=10e-9,
                              dXYZ_factors=None, fXYZ_factors=None,
                              reduce_list=False, return_calibration=False,
                              return_settings=False):
        return self.data_pairs(i=i, axis=axis, direction=direction,
                               decimate=decimate, fe_data=True, twoD=twoD,
                               posmin=posmin, dXYZ_factors=dXYZ_factors,
                               fXYZ_factors=fXYZ_factors,
                               reduce_list=reduce_list,
                               return_calibration=return_calibration,
                               return_settings=return_settings)

    def data_pairs(self, i=None, axis=None, direction=None, decimate=None,
                   fe_data=False, twoD=False, posmin=10e-9, dXYZ_factors=None,
                   fXYZ_factors=None, reduce_list=False,
                   return_calibration=False, return_settings=False):
        """
        Return a dictionary of force extension values of stress release pairs.

        Parameters
        ----------
        i : int
            Index of the force extension pair to be yielded.
        axis : str, optional
            See method `self.stress_release_pairs()`.
        direction : str, optional
            See method `self.stress_release_pairs()`.
        decimate : int, optional
            See method `seld.stress_release_pairs()`.
        posmin : float
            The `posmin` is used to decide wether the magnitude of the force
            has to be corrected with the sign depending on the position. The
            `posmin` sets the value the position signal has to be deflected at
            least to be counted as active pulling on the bead. The value should
            at least be >= 12 times the standard deviation of the unexcited
            position signal.
            Smaller values could (depending on the number of datapoints)
            possibly lead to falsly detected excitation of the signal.

        Returns
        ------
        dict
            1D numpy.ndarray of type float
                Extension values of stress cycles in m.
            1D numpy.ndarray of type float
                Force values of stress cycles in N.
            1D numpy.ndarray of type str
                (str, str, str), containing the axis, direction, and the cycle.
                Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
                Cycle is 'stress'.
            1D numpy.ndarray of type float
                Extension values of release cycles in m.
            1D numpy.ndarray of type float
                Force values of release cycles in N.
            1D numpy.ndarray of type str
                (str, str, str), containing the axis, direction, and the cycle.
                Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
                Cycle is 'release'.
        """
        pairs = self.stress_release_pairs(i=i, axis=axis, direction=direction,
                                          decimate=decimate, slices=True,
                                          reduce_list=False)

        # Get time, extension, and force for all stress/release pairs. Set the
        # start to the first stress and the stop to the last release cycle
        start = pairs['stress']['idx'][0].start
        stop = pairs['release']['idx'][-1].stop
        samples = slice(start, stop)

        t = self.timevector[samples]
        keys = ['psdXYZ', 'positionXYZ']
        if fe_data:
            _data = self.force_extension(samples=samples, twoD=twoD,
                                         posmin=posmin,
                                         dXYZ_factors=dXYZ_factors,
                                         fXYZ_factors=fXYZ_factors,
                                         return_calibration=return_calibration)
            keys = keys + ['extension', 'force']
        else:
            _data = self.raw_data(samples=samples,
                                  return_calibration=return_calibration)

        # Separate all data into stress/release pairs
        data = {}

        for cycle in ['stress', 'release']:
            data[cycle] = {}
            data[cycle]['time'] = []
            for key in keys:
                data[cycle][key] = []
            idcs = pairs[cycle]['idx']
            for idx in idcs:
                i = slice(idx.start - start, idx.stop - start, idx.step)
                data[cycle]['time'].append(t[i])
                for key in keys:
                    data[cycle][key].append(_data[key][i])

        # If only one pair, remove list and select only pair
        if reduce_list and len(data['stress']['time']) == 1:
            for cycle in ['stress', 'release']:
                for key in data[cycle]:
                    data[cycle][key] = data[cycle][key][0]

        if return_calibration:
            data['calibration'] = _data['calibration']
        if return_settings:
            data['settings'] = {
                'i': i,
                'axis': axis,
                'direction': direction,
                'decimate': decimate,
                'twoD': twoD,
                'posmin': posmin,
                'dXYZ_factors': dXYZ_factors,
                'fXYZ_factors': fXYZ_factors
            }

        return data


    def samples(self, i=None, cycle=None, axis=None, direction=None,
                decimate=None):
        """
        Get index of samples to be used for functions displacement, force, etc.

        Parameters
        ----------
        cycle : str or list of str
            cycle can be either 'stress' or 'release'
        """
        cycle = ['stress', 'release'] if cycle is None else cycle
        start = 'stress'
        stop = 'release'
        if 'stress' in cycle and not 'release' in cycle:
            # stress only
            stop = start
        elif 'release' in cycle and not 'stress' in cycle:
            # release only
            start = stop

        pairs = self.stress_release_pairs(i=i, axis=axis, direction=direction,
                                          decimate=decimate, slices=True)

        idx_start = pairs[start]['idx'][0].start
        idx_stop = pairs[stop]['idx'][-1].stop
        idx = slice(idx_start, idx_stop, decimate)
        return idx


    def get_data(self, i=None, cycle=None, axis=None, direction=None,
                 decimate=None, samples=None, **kwargs):
        if samples is None and not (i is None and cycle is None
                                    and axis is None and direction is None
                                    and decimate is None):
            samples = self.samples(i=i, cycle=cycle, axis=axis,
                                   direction=direction, decimate=decimate)
        data = super().get_data(**kwargs, samples=samples)
        return data


    def displacementXYZ(self, samples=None, dXYZ_factors=None):
        """
        Displacement in m with height dependent calibration factors for X, Y
        and Z.
        """
        # Get extension (in a fast way)
        data = self.get_data(traces=['psdXYZ', 'positionZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionZ = data[:, [3]]
        calibration = self.calibration

        dispXYZ = displacementXYZ(psdXYZ, calibration=calibration,
                                  positionZ=positionZ,
                                  dXYZ_factors=dXYZ_factors)
        return dispXYZ


    def forceXYZ(self, samples=None, dXYZ_factors=None, fXYZ_factors=None):
        """
        Force in N, that is acting on the tether
        """
        data = self.get_data(traces=['psdXYZ', 'positionZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionZ = data[:, [3]]
        calibration = self.calibration

        fXYZ = forceXYZ(psdXYZ, calibration=calibration, positionZ=positionZ,
                        dXYZ_factors=dXYZ_factors, fXYZ_factors=fXYZ_factors)
        return fXYZ


    def force(self, samples=None, twoD=False, posmin=10e-9, dXYZ_factors=None,
              fXYZ_factors=None):
        """
        Magnitude of the force in N acting on the tethered molecule (1D
        numpy.ndarray).

        Parameters
        ----------
        posmin : float
            The `posmin` is used to decide wether the magnitude of the force
            has to be corrected with the sign depending on the position. The
            `posmin` sets the value the position signal has to be deflected at
            least to be counted as active pulling on the bead. The value should
            at least be >= 12 times the standard deviation of the unexcited
            position signal.
            Smaller values could (depending on the number of datapoints)
            possibly lead to falsly detected excitation of the signal.
        """
        # Get force (in a fast way)
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionXY = data[:, 3:5]
        positionZ = data[:, [5]]
        calibration = self.calibration

        # 2D or 3D calculation of the distance in Z
        if twoD:
            psdXYZ[:, Z] = 0.0

        fXYZ = forceXYZ(psdXYZ, calibration=calibration, positionZ=positionZ,
                        dXYZ_factors=dXYZ_factors, fXYZ_factors=fXYZ_factors)
        f = force(fXYZ, positionXY, posmin=posmin)
        return f


    def distanceXYZ(self, samples=None, dXYZ_factors=None):
        """
        Distance of the attachment point to the bead center as a 3D vector.
        """
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionXYZ = data[:, 3:6]
        calibration = self.calibration

        distXYZ = distanceXYZ(positionXYZ, psdXYZ=psdXYZ,
                              calibration=calibration,
                              dXYZ_factors=dXYZ_factors)
        return distXYZ


    def distance(self, samples=None, twoD=False, posmin=10e-9,
                 dXYZ_factors=None):
        """
        Distance of the attachment point to the bead center.

        Parameters
        ----------
        posmin : float
            The `posmin` is used to decide wether the magnitude of the force
            has to be corrected with the sign depending on the position. The
            `posmin` sets the value the position signal has to be deflected at
            least to be counted as active pulling on the bead. The value should
            at least be >= 12 times the standard deviation of the unexcited
            position signal.
            Smaller values could (depending on the number of datapoints)
            possibly lead to falsly detected excitation of the signal.
        """
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionXYZ = data[:, 3:6]
        positionXY = data[:, 3:5]
        calibration = self.calibration

        # 2D or 3D calculation of the distance in Z
        if twoD:
            psdXYZ[:, Z] = 0.0

        distXYZ = distanceXYZ(positionXYZ, psdXYZ=psdXYZ,
                              calibration=calibration,
                              dXYZ_factors=dXYZ_factors)
        dist = distance(distXYZ, positionXY, posmin=posmin)
        return dist


    def extension(self, samples=None, twoD=False, posmin=10e-9):
        """
        Extension of the tethered molecule in m.

        Parameters
        ----------
        posmin : float
            The `posmin` is used to decide wether the magnitude of the force
            has to be corrected with the sign depending on the position. The
            `posmin` sets the value the position signal has to be deflected at
            least to be counted as active pulling on the bead. The value should
            at least be >= 12 times the standard deviation of the unexcited
            position signal.
            Smaller values could (depending on the number of datapoints)
            possibly lead to falsly detected excitation of the signal.
        """
        calibration = self.calibration

        dist = self.distance(samples=samples, twoD=twoD, posmin=posmin)
        e = extension(dist, calibration.radius)
        return e


    def raw_data(self, i=None, cycle=None, axis=None, direction=None,
                 samples=None, return_calibration=False,
                 return_settings=False):
        # Get extension and force (in a fast way)
        data = self.get_data(i=i, cycle=cycle, axis=axis, direction=direction,
                             traces=['psdXYZ', 'positionXYZ'], samples=samples)
        psdXYZ = data[:, 0:3]
        positionXYZ = data[:, 3:6]

        raw_data = {
            'psdXYZ': psdXYZ,
            'positionXYZ': positionXYZ
        }
        if return_calibration:
            calibration = self.calibration
            raw_data['calibration'] = {
                'beta0': calibration._beta[0],
                'betam': calibration._beta[1],
                'kappa0': calibration._kappa[0],
                'kappam': calibration._kappa[1],
                'dsurf': calibration.dsurf,
                'psurf': calibration.psurf,
                'touchdown': calibration.touchdown,
                'radius': calibration.radius,
                'correct_radius': calibration.correct_radius,
                'corrfactor': calibration.calibsource.corrfactor,
                'radiusspec': calibration.calibsource.radiusspec,
                'focalshift': calibration.focalshift
            }
        if return_settings:
            raw_data['settings'] = {
                'i': i,
                'cycle': cycle,
                'axis': axis,
                'direction': direction,
                'samples': samples,
            }
        return raw_data


    def force_extension(self, i=None, cycle=None, axis=None, direction=None,
                        samples=None, twoD=False, posmin=10e-9,
                        dXYZ_factors=None, fXYZ_factors=None,
                        return_calibration=False, return_settings=False):
        """
        Extension (m, first column) of and force (N, second column) acting
        on the tethered molecule (2D numpy.ndarray).

        Parameters
        ----------
        posmin : float
            The `posmin` is used to decide wether the magnitude of the force
            has to be corrected with the sign depending on the position. The
            `posmin` sets the value the position signal has to be deflected at
            least to be counted as active pulling on the bead. The value should
            at least be >= 12 times the standard deviation of the unexcited
            position signal.
            Smaller values could (depending on the number of datapoints)
            possibly lead to falsly detected excitation of the signal.
        """
        data = self.raw_data(i=i, cycle=cycle, axis=axis,
                             direction=direction, samples=samples,
                             return_calibration=return_calibration)
        psdXYZ = data['psdXYZ']
        positionXYZ = data['positionXYZ']
        positionXY = positionXYZ[:, 0:2]
        positionZ = positionXYZ[:, [2]]
        calibration = self.calibration

        # 2D or 3D calculation of the distance in Z
        if twoD:
            psdXYZ[:, Z] = 0.0
        dispXYZ = displacementXYZ(psdXYZ, calibration=calibration,
                                  positionZ=positionZ,
                                  dXYZ_factors=dXYZ_factors)
        distXYZ = distanceXYZ(positionXYZ, dispXYZ=dispXYZ,
                              calibration=calibration,
                              dXYZ_factors=dXYZ_factors)
        dist = distance(distXYZ, positionXY, posmin=posmin)
        e = extension(dist, calibration.radius)

        fXYZ = forceXYZ(psdXYZ, calibration=calibration, positionZ=positionZ,
                        dXYZ_factors=dXYZ_factors, fXYZ_factors=fXYZ_factors)
        f = force(fXYZ, positionXY, posmin=posmin)

        data.update({
            'extension': e,
            'force': f
        })
        if return_settings:
            data['settings'] = {
                'i': i,
                'cycle': cycle,
                'axis': axis,
                'direction': direction,
                'samples': samples,
                'twoD': twoD,
                'posmin': posmin,
                'dXYZ_factors': dXYZ_factors,
                'fXYZ_factors': fXYZ_factors
            }
        return data


    def info(self, i=0):
        print_info(self, i=i)

    @property
    def angle(self):
        """
        Returns a dictionary of angles, for three corners (A), (B), and (C),
        calculated by (F)orce and (D)istance:
            B
            |\
          a | \ c
            |  \
            |___\
           C  b  A
        """
        return angle(self.force, self.forceXYZ, self.excited_axis,
                     self.distance, self.distanceXYZ)

    @property
    def rightstress(self):
        return self.sections(direction='right', cycle='stress',
                             range_concat=True)

    @property
    def leftstress(self):
        return self.sections(direction='left', cycle='stress',
                             range_concat=True)

    @property
    def leftrelease(self):
        return self.sections(direction='left', cycle='release',
                             range_concat=True)

    @property
    def rightrelease(self):
        return self.sections(direction='right', cycle='release',
                             range_concat=True)

    @property
    def stress(self):
        return self.sections(cycle='stress', range_concat=True)

    @property
    def release(self):
        return self.sections(cycle='release', range_concat=True)

    @property
    def rightstresses(self):
        return self.sections(direction='right', cycle='stress')

    @property
    def leftstresses(self):
        return self.sections(direction='left', cycle='stress')

    @property
    def leftreleases(self):
        return self.sections(direction='left', cycle='release')

    @property
    def rightreleases(self):
        return self.sections(direction='right', cycle='release')

    @property
    def stresses(self):
        return self.sections(cycle='stress')

    @property
    def releases(self):
        return self.sections(cycle='release')


# Define constants for convenient handling
X = 0
Y = 1
Z = 2
XY = hp.slicify([X, Y])
XZ = hp.slicify([X, Z])
YZ = hp.slicify([Y, Z])
XYZ = hp.slicify([X, Y, Z])


def print_info(tether, i=0):
    pairs = tether.stress_release_pairs(i=i, slices=True)
    stress = pairs['stress']['idx'][0]
    release = pairs['release']['idx'][0]
    stress_info = pairs['stress']['info'][0]
    release_info = pairs['release']['info'][0]
    resolution = tether.resolution
    start_stress_t = stress.start / resolution
    stop_stress_t = stress.stop / resolution
    start_release_t = release.start / resolution
    stop_release_t = release.stop / resolution
    print("Stress release pair: #{:03d}".format(i))
    print("Focal shift: {:.3f}".format(tether.calibration.focalshift))
    print("Axis: {}".format(stress_info[0][0]))
    print("  stress")
    print("    t:  {:.3f} s".format(stop_stress_t - start_stress_t))
    print("        {:.2f} s - {:.2f} s".format(start_stress_t, stop_stress_t))
    print("    datapoints: {}".format(stress.stop - stress.start))
    print("    z0: {:.2f} nm".format(
        - np.median(tether.get_data('positionZ', stress)) * 1e9))
    print("    h0: {:.2f} nm".format(
        - np.median(tether.get_data('positionZ', stress)) * 1e9
                                    * tether.calibration.focalshift))
    print("  release")
    print("    t:  {:.3f} s".format(stop_release_t - start_release_t))
    print("        {:.2f} s - {:.2f} s".format(start_release_t,
                                               stop_release_t))
    print("    datapoints: {}".format(release.stop - release.start))
    print("    z0: {:.2f} nm".format(
        - np.median(tether.get_data('positionZ', release)) * 1e9))
    print("    h0: {:.2f} nm".format(
        - np.median(tether.get_data('positionZ', release)) * 1e9
                                    * tether.calibration.focalshift))


def displacementXYZ(psdXYZ, calibration=None, positionZ=0.0, beta=None,
                    dXYZ_factors=None):
    """
    Displacement in m with height dependent calibration factors for X, Y
    and Z.
    """
    if calibration is not None:
        dispXYZ = calibration.displacement(psdXYZ, positionZ=positionZ)
    else:
        dispXYZ = psdXYZ * beta

    # Optionally use correctur factors for the calculated displacement
    dXYZ_factors = 1 if dXYZ_factors is None else dXYZ_factors

    dispXYZ *= dXYZ_factors

    return dispXYZ


def forceXYZ(psdXYZ=None, dispXYZ=None, calibration=None, positionZ=0.0,
             beta=None, kappa=None, dXYZ_factors=None, fXYZ_factors=None):
    """
    Force acting on the bead
    """
    if dispXYZ is None:
        dispXYZ = displacementXYZ(psdXYZ, calibration=calibration,
                                  positionZ=positionZ, beta=beta,
                                  dXYZ_factors=dXYZ_factors)

    if calibration is not None:
        fXYZ = calibration.force(dispXYZ, positionZ=positionZ)
    else:
        fXYZ = dispXYZ * kappa

    # Optionally use correctur factors for the calculated force
    fXYZ_factors = 1 if fXYZ_factors is None else fXYZ_factors

    # A displacement of the bead in the positive (negative) direction results
    # in a force acting on the bead in the opposite negative (positive)
    # direction.
    fXYZ *= - 1.0 * fXYZ_factors

    return fXYZ


def force(forceXYZ, positionXY, posmin=10e-9, sign_axes=None):
    """
    Parameters
    ----------
    forceXYZ : 2D numpy.ndarray of type float
        forceXYZ.shape[1] can consist of either 3 (XYZ) or 2 (XY) axes
    posmin : float
        The `posmin` is used to decide wether the magnitude of the force has to
        be corrected with the sign depending on the position. The `posmin` sets
        the value the position signal has to be deflected at least to be
        counted as active pulling on the bead. The value should at least be >=
        12 times the standard deviation of the unexcited position signal.
        Smaller values could (depending on the number of datapoints) possibly
        lead to falsly detected excitation of the signal.
    """
    if posmin is None:
        posmin = 0
    # The sign of the magnitude, i.e. the direction of the force is important
    # for the noise around +/- 0 N.
    # The sign of the magnitude of the forceXY depends on the positionXY. If we
    # pull the bead to one side and the bead is displaced to the same side, we
    # get a positive magnitude of the force. If the bead is displaced to the
    # opposite side of the one we are pulling to, we get a negative magnitude
    # of the force. Keep in mind that a positive displacement results in an
    # opposite directed negative force acting on the bead and vice versa. The
    # sign of the magnitude of the forceZ is independent of the positionZ.
    signF = np.sign(forceXYZ)
    # The position determines the direction of the force only of the axes where
    # we actively pull on the bead
    posmax = np.max(np.abs(positionXY), axis=0)
    idx = posmax >= posmin
    signF[:, XY][:, idx] *= - np.sign(positionXY[:, idx])

    # Square the forces and account for the signs
    force_sq = forceXYZ**2
    force_sq_sum = np.sum(force_sq, axis=1)
    # Calculate a "weighted" sign. Greater forces have greater influence on the
    # final sign
    axs = XY if sign_axes is None else sign_axes
    sign_force_sum = np.sign(np.sum(force_sq[:,axs] * signF[:,axs], axis=1))
    force = np.sqrt(force_sq_sum) * sign_force_sum
    return force


def distanceXYZ(positionXYZ, psdXYZ=None, dispXYZ=None, calibration=None,
                beta=None, radius=None, focalshift=None, clip_Z=True,
                dXYZ_factors=None):
    """
    Distance of the attachment point to the bead center as a 3D vector.
    positionXYZ, displacementXYZ and radius need to have the same unit.
    If radius is 0.0, the positionZ is not corrected by the radius. You should
    set radius to 0.0, if the positionZ is defined in such a way, that
    positionZ is 0.0, where the bead center would be on the glass surface.
    If focalshift is 1.0, the positionXYZ is not corrected by the focalshift.

    Parameters
    ----------
    positionXYZ : 2D np.array of type float
    displacementXYZ : 2D np.array of type float
    radius : float
    focalshift : float
    clip_Z : bool
        The distance of the attachment point to the center of the bead cannot
        be smaller than the radius. Therefore, clip the data to be at least as
        great as the radius.
        However, values much smaller than the radius could indicate an errornes
        calibration with too small displacement sensitivities, which would lead
        to too small displacements in Z and in turn to negative distances.
        Therefore, if you want to check for this kind of error, switch off the
        clip_Z functionality.
    """
    if focalshift is None:
        focalshift = calibration.focalshift
    if radius is None:
        radius = calibration.radius
    positionZ = positionXYZ[:,2]

    if dispXYZ is None:
        dispXYZ = displacementXYZ(psdXYZ, calibration=calibration,
                                  positionZ=positionZ, beta=beta,
                                  dXYZ_factors=dXYZ_factors)

    # distance, point of attachment of DNA
    # displacement, displacement of bead out of trap center radius
    distanceXYZ = positionXYZ.copy()
    # distance from attachment point to center of bead
    # attachmentXY - displacementXY
    distanceXYZ[:, 0:2] -= dispXYZ[:, 0:2]
    # Direction of distance vector is from attachment point to the bead center
    distanceXYZ[:, 0:2] *= -1

    # If the bead is free (i.e. above the surface with a distance Z > 0), a
    # movement of the positionZ leads to a distance change reduced by the focal
    # shift.
    # If the bead is on the surface, a movement of the positionZ leads to a
    # distance change independant of the focal shift.
    idx_free = positionXYZ[:, 2] < 0
    idx_touch = positionXYZ[:, 2] >= 0
    distanceXYZ[idx_free, 2] = (- positionXYZ[idx_free, 2] * focalshift
                                + radius
                                # distanceZ + radius + displacementZ
                                + dispXYZ[idx_free, 2])
    distanceXYZ[idx_touch, 2] = (- positionXYZ[idx_touch, 2]
                                 + radius
                                 # distanceZ + radius + displacementZ
                                 + dispXYZ[idx_touch, 2])

    if clip_Z:
        distanceXYZ[:, 2] = distanceXYZ[:, 2].clip(min=radius)
    # A positive positionZ signal (positionZ upwards) corresponds to a
    # decreasing (negative) distance of the bead to the surface:
    #   -> distanceZ ~ - positionZ
    # A movement of the positionZ gets reduced by the focalshift:
    #   -> distanceZ = - positionZ * focalshift
    # The distanceZ (positionZ) is 0, where the bead touches the surface:
    #   -> center of the bead is at distanceZ + radius
    # The bead is stressed down with increasing positive distanceZ. This leads
    # to a negative displacement, which reduces the distance:
    #   -> distanceZ + displacement

    # distance from attachment point to bead center
    return distanceXYZ


def displacement(displacementXYZ, positionXY, posmin=10e-9):
    """
    Calculate the displacemnet of the microsphere.

    Parameters
    ----------
    displacementXYZ : 2D numpy.ndarray of type float
        displacementXYZ.shape[1] can consist of either 3 (XYZ) or 2 (XY) axes
    posmin : float
        The `posmin` is used to decide wether the magnitude of the force has to
        be corrected with the sign depending on the position. The `posmin` sets
        the value the position signal has to be deflected at least to be
        counted as active pulling on the bead. The value should at least be >=
        12 times the standard deviation of the unexcited position signal.
        Smaller values could (depending on the number of datapoints) possibly
        lead to falsly detected excitation of the signal.
    """
    return - distance(displacementXYZ, positionXY, posmin=posmin)


def distance(distanceXYZ, positionXY, posmin=10e-9, sign_axes=None):
    """
    Calculate the distance of the attachment point to the bead center.

    Parameters
    ----------
    distanceXYZ : 2D numpy.ndarray of type float
        distanceXYZ.shape[1] can consist of either 3 (XYZ) or 2 (XY) axes
    posmin : float
        The `posmin` is used to decide wether the magnitude of the force has to
        be corrected with the sign depending on the position. The `posmin` sets
        the value the position signal has to be deflected at least to be
        counted as active pulling on the bead. The value should at least be >=
        12 times the standard deviation of the unexcited position signal.
        Smaller values could (depending on the number of datapoints) possibly
        lead to falsly detected excitation of the signal.
    """
    if posmin is None:
        posmin = 0
    # The sign of the magnitude, i.e. the direction of the distance is
    # important for the noise around +/- 0 N.
    # The sign of the magnitude of the distanceXY depends on the positionXY. If
    # we pull the bead to one side and the bead is displaced to the same side,
    # we get a positive magnitude of the distance. If the bead is displaced to
    # the opposite side of the one we are pulling to, we get a negative
    # magnitude of the distance. Keep in mind that a positive displacement
    # results in an opposite directed pointing distance. The sign of the
    # magnitude of the distanceZ is independent of the positionZ.
    signD = np.sign(distanceXYZ)
    # The position determines the direction of the force only of the axes where
    # we actively pull on the bead
    posmax = np.max(np.abs(positionXY), axis=0)
    idx = posmax >= posmin
    signD[:, XY][:, idx] *= - np.sign(positionXY[:, idx])

    # Square the distances and account for the signs
    distance_sq = distanceXYZ**2
    dist_sq_sum = np.sum(distance_sq, axis=1)
    # Calculate a "weighted" sign. Greater distances have greater influence on
    # the final sign
    axs = XYZ if sign_axes is None else sign_axes
    sign_dist_sum = np.sign(np.sum(distance_sq[:,axs] * signD[:,axs], axis=1))
    return np.sqrt(dist_sq_sum) * sign_dist_sum


def extension(distance, radius):
    """
    Calculate the extension of the DNA by simply subtracting the radius from
    the distance.
    """
    return distance - radius


def angle(force, forceXYZ, excited_axis, distance, distanceXYZ):
    """
    Returns a dictionary of angles, for three corners (A), (B), and (C),
    calculated by (F)orce and (D)istance:
        B
        |\
      a | \ c
        |  \
        |___\
       C  b  A
    """

    Fabs = force
    Fz = forceXYZ[:, 2]
    if excited_axis == 'X':
        ea = 0
    else:
        ea = 1
    Fxy = forceXYZ[:, ea]

    dist = distance
    dist3D = distanceXYZ
    Dz = dist3D[:, 2]
    Dxy = dist3D[:, ea]

    AF = _angle(Fz, Fxy, Fabs, angle='A')
    BF = _angle(Fz, Fxy, Fabs, angle='B')
    CF = _angle(Fz, Fxy, Fabs, angle='C')
    AD = _angle(Dz, Dxy, dist, angle='A')
    BD = _angle(Dz, Dxy, dist, angle='B')
    CD = _angle(Dz, Dxy, dist, angle='C')

    angle = dict(list(zip(['AF', 'BF', 'CF', 'AD', 'BD', 'CD'],
                          [AF, BF, CF, AD, BD, CD])))

    return angle


def _angle(a, b, c, angle='A'):
    """
        B
        |\
      a | \ c
        |  \
        |___\
       C  b  A
    """
    if angle == 'A':
        _a = b
        _b = c
        _c = a
    elif angle == 'B':
        _a = a
        _b = c
        _c = b
    else:
        _a = a
        _b = b
        _c = c

    cos_angle = (_a**2 + _b**2 - _c**2) / (2 * _a * _b)
    return np.arccos(np.fabs(cos_angle)) * 180.0 / np.pi
