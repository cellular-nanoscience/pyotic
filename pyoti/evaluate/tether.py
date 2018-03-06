# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 13:36:39 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import interact, IntSlider

from . import signal as sn
from . import dna
from . import evaluate
from .. import gui
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
    def __init__(self, region=None, calibration=None, resolution=1000,
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

        traces_sf = kwargs.pop('traces_sf', 'positionXY')

        # resolution for (almost) all properties, e.g. displacement, force, ...
        # resolution 1.000 Hz, samplingrate 40.000 Hz -> 40 points
        super().__init__(region=region, calibration=calibration,
                         resolution=resolution, filter_time=filter_time,
                         sf_class=CycleSectioner, traces_sf=traces_sf,
                         resolution_sf=resolution_sf,
                         filter_time_sf=filter_time_sf)

        self.fe_figure = None

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

    def stress_release_pairs(self, axis=None, direction=None, i=None,
                             slices=True, decimate=None, info=False, **kwargs):
        """
        Calculate start/stop indices (as segments or slices) of stress/release
        cycle pairs of the requested axis and direction.

        Parameters
        ----------
        axis : str or list of str, optional
            'x', 'y', or ['x', 'y'] (default)
        direction : str or list of str, optional
            'left', 'right', or ['left', 'right'] (default)
        i : int, optional
            The index of the stress release pair to be selected and returned
            from all stress release pairs calculated. Depending on `axis` and
            `direction` the total availabe number of stress release pairs may
            vary.
        slices : bool, optional
            Set to False to return segments instead of slices. Default is True.
        decimate : int, optional
            Used to set the step attribute of the returned slices. Only
            evaluated, if `slices` is True.
        info : bool, optional
        **kwargs : dict, optional
            Only used for compatibility purposes, to be able to call method
            with parameters not defined in method definition.

        Returns
        -------
        stresses : 1D np.ndarray of slices or segments
        releases : 1D np.ndarray of slices or segments
        stress_infos : 1D np.ndarray of 1D np.ndarrays of type str
            The str arrays have the form of [axis, direction, cycle]. Only
            returned, if parameter `info` is True.
        release_infos : 1D np.ndarray of 1D np.ndarrays of type str
            The str arrays have the form of [axis, direction, cycle]. Only
            returned, if parameter `info` is True.
        """
        if axis is None:
            axis = ['x', 'y']
        if direction is None:
            direction = ['left', 'right']
        if isinstance(direction, str):
            direction = [direction]
        stress_segments, _stress_infos = self.sections(axis=axis,
                                                       direction=direction,
                                                       cycle='stress',
                                                       slices=False,
                                                       info=True)
        release_segments, _release_infos = self.sections(axis=axis,
                                                         direction=direction,
                                                         cycle='release',
                                                         slices=False,
                                                         info=True)
        extrema = self._extrema(axis=axis)
        stresses = np.empty((0, 2), dtype=int)
        releases = np.empty((0, 2), dtype=int)
        stress_infos = np.empty((0, 3), dtype=str)
        release_infos = np.empty((0, 3), dtype=str)

        # Group all stress/release cycle pairs according to the extrema
        # A stress release pair corresponds to one extremum, only if the stop
        # of the stress and the start of the release segment equal the
        # extremum. However, either one of the stress or release segment can be
        # missing.
        for ext in extrema:
            # Find stress whose stop is equal to extremum
            stress_idx = stress_segments[:, 1] == ext
            stress = stress_segments[stress_idx]
            stress_info = _stress_infos[stress_idx]
            # Find release whose start is equal to extremum
            release_idx = release_segments[:, 0] == ext
            release = release_segments[release_idx]
            release_info = _release_infos[release_idx]
            # Create stress slice, if stress is empty, but release is valid
            if stress.size == 0 and release.size > 0:
                stress = np.array([[ext, ext]])
            # Create release slice, if release is empty, but stress is valid
            if release.size == 0 and stress.size > 0:
                release = np.array([[ext, ext]])
            if stress.size > 0 and release.size > 0:
                # p_p = np.r_[stress, release]  # concatenate stress/release
                # cycle pair
                # concatenate with other cycles
                stresses = np.r_[stresses, stress]
                releases = np.r_[releases, release]
                stress_infos = np.r_[stress_infos, stress_info]
                release_infos = np.r_[release_infos, release_info]

        # Convert segments into slices
        if slices:
            stresses = sn.idx_segments_to_slices(stresses, decimate=decimate)
            releases = sn.idx_segments_to_slices(releases, decimate=decimate)

        # Get the maximum number of stress release pairs
        stop = len(stresses)
        # Prevent index overflow and allow negative indices
        if i is not None:
            if i < 0:
                i = stop + i
            i = max(0, i)
            i = min(i, stop - 1)
            s = slice(i, i + 1)
        else:
            s = slice(0, stop)

        if info:
            return stresses[s], releases[s], stress_infos[s], release_infos[s]
        else:
            return stresses[s], releases[s]

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
            stresses, releases = self.stress_release_pairs(axis=axis,
                                                           slices=False)
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
            s = self.get_data(trace)
            rstr, rrls = self.stress_release_pairs(axis=axis,
                                                   direction='right')
            lstr, lrls = self.stress_release_pairs(axis=axis, direction='left')

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

        ax.set_ylabel("Signal positionX and Y (um)")
        ax.set_xlabel("Time (s)")
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

    def force_extension_pair(self, i, axis=None, direction=None, decimate=None,
                             time=False, twoD=False):
        """
        Calculate the force extension pair with index `i`.

        Parameters
        ----------
        i : int
            Index of the force_extension pair to be returned.
        axis : str
            See method `self.stress_release_pairs()`.
        direction : str
            See method `self.stress_release_pairs()`.
        decimate : int
            See method `self.stress_release_pairs()`.

        Returns
        -------
        1D numpy.ndarray of type float
            Extension values of stress cycles in nm.
        1D numpy.ndarray of type float
            Force values of stress cycles in pN.
        1D numpy.ndarray of type str
            (str, str, str), containing the axis, direction, and the cycle.
            Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
            Cycle is 'stress'.
        1D numpy.ndarray of type float
            Extension values of release cycles in nm.
        1D numpy.ndarray of type float
            Force values of release cycles in pN.
        1D numpy.ndarray of type str
            (str, str, str), containing the axis, direction, and the cycle.
            Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
            Cycle is 'release'.
        """
        fe_pair = self.force_extension_pairs(axis=axis, direction=direction,
                                             i=i, decimate=decimate, time=time,
                                             twoD=twoD)
        return next(fe_pair)

    def force_extension_pairs(self, axis=None, direction=None, i=None,
                              decimate=None, time=False, twoD=False):
        """
        Return a generator for force extension values of stress release pairs.

        Parameters
        ----------
        str_rls_pairs : tuple of 4 numpy.ndarrays, optional
        axis : str, optional
            See method `self.stress_release_pairs()`.
        direction : str, optional
            See method `self.stress_release_pairs()`.
        i : int
            Index of the force extension pair to be yielded.
        decimate : int, optional
            See method `seld.stress_release_pairs()`.

        Yields
        ------
        1D numpy.ndarray of type float
            Extension values of stress cycles in nm.
        1D numpy.ndarray of type float
            Force values of stress cycles in pN.
        1D numpy.ndarray of type str
            (str, str, str), containing the axis, direction, and the cycle.
            Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
            Cycle is 'stress'.
        1D numpy.ndarray of type float
            Extension values of release cycles in nm.
        1D numpy.ndarray of type float
            Force values of release cycles in pN.
        1D numpy.ndarray of type str
            (str, str, str), containing the axis, direction, and the cycle.
            Axis can be either 'x' or 'y'. Direction can be 'left', or 'right'.
            Cycle is 'release'.
        """
        str_rls_pairs = self.stress_release_pairs(axis=axis,
                                                  direction=direction,
                                                  i=i,
                                                  decimate=decimate,
                                                  slices=True,
                                                  info=True)
        strs, rlss, stri, rlsi = str_rls_pairs

        # Get the start from the first stress and the stop from the last
        # release cycle
        start = strs[0].start
        stop = rlss[-1].stop
        samples = slice(start, stop)

        # Get extension, force, and stress/release pairs
        e_f = self._force_extension(samples=samples, twoD=twoD) * 1000  # nm,pN
        e = e_f[:, 0]
        f = e_f[:, 1]

        # Get the time
        if time:
            t = self.timevector[samples]

        # Yield all stress/release extension/force data pairs
        for st, rl, sti, rli in zip(strs, rlss, stri, rlsi):
            _st = slice(st.start - start, st.stop - start, st.step)
            _rl = slice(rl.start - start, rl.stop - start, rl.step)
            est = e[_st]
            fst = f[_st]
            erl = e[_rl]
            frl = f[_rl]
            if time:
                tst = t[_st]
                trl = t[_rl]
                yield est, fst, sti, erl, frl, rli, tst, trl
            else:
                yield est, fst, sti, erl, frl, rli

    def init_fe_fig(self, show=True, autolimit=True, plot_params=None):
        """
        Initialize the figure where the force extension pairs should be plotted
        to.

        Parameters
        ----------
        show : bool, optional
            Show the figure in the notebook after it has been initialized.
        autolimit : bool, optional
            Set the limits of the force extension plot according to the min/max
            values of the force and extension.
        plot_params : dict, optional
        """
        xlim = None
        ylim = None
        if self.fe_figure is not None:
            ax = self.fe_figure.axes[0]
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            plt.close(self.fe_figure)

        # Initialize proper plot parameters
        gui.set_plot_params(plot_params=plot_params)

        # Create static figure
        self.fe_figure = _force_extension_figure()

        if autolimit:
            self.autolimit()
        else:
            ax = self.fe_figure.axes[0]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        if show:
            self.fe_figure.show()

    def autolimit(self, samples=None, e=None, f=None, xlim=None, ylim=None,
                  set_limits=True):
        """
        Determine xlim and ylim values for the force extension pair plots.

        Parameters
        ----------
        samples : int, slice or index array, optional
            Samples to get extension and force from.
        e : 1D numpy.ndarray of floats, optional
            Extension in nm. Takes precedence over extension determined with
            `samples`.
        f : 1D numpy.ndarray of floats, optional
            Force in pN. Takes precedence over force determined with `samples`.
        xlim : (float, float), optional
            Xlimit of force extension axis. Takes precedence over xlim
            determined with `e`.
        ylim : (float, float), optional
            Ylimit of force extension axis. Takes precedence over ylim
            determined with `f`.
        set_limits : bool, optional
            Set xlim and ylim of the force extension plot.

        Returns
        -------
        (float, float)
            The xlim
        (float, float)
            The ylim
        """
        if samples is None \
                and (xlim is None and e is None) \
                or (ylim is None and f is None):
            # Get the start/stop indices of the data to be used to determine
            # the min max values
            sts, rls = self.stress_release_pairs()
            start = sts[0].start
            stop = rls[-1].stop
            samples = slice(start, stop)

        if xlim is None and ylim is None and e is None and f is None:
                e_f = self._force_extension(samples=samples) * 1000  # nm, pN
                e = e_f[:, 0]
                f = e_f[:, 1]
        if xlim is None and e is None:
            e = self._extension(samples=samples) * 1000  # nm
        if ylim is None and f is None:
            f = self._force(samples=samples) * 1000  # pN

        if xlim is None:
            e_min = e.min()
            e_max = e.max()
            e_diff = (e_max - e_min) * 0.02
            xlim = (e_min - e_diff, e_max + e_diff)

        if ylim is None:
            f_min = f.min()
            f_max = f.max()
            f_diff = (f_max - f_min) * 0.02
            ylim = (f_min - f_diff, f_max + f_diff)

        # Set the limits
        if set_limits:
            fig = self.fe_figure
            ax = fig.axes[0]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            fig.canvas.draw()

        # Return the set limits
        return xlim, ylim

    def update_force_extension_plot(self, fe_pair=None, xlim=None, ylim=None,
                                    title=None, draw=True, bps=None):
        """
        Update the figure with force extension data.

        Parameters
        ----------
        fe_pair : tuple of 6 1D numpy.ndarrays, optional
            As returned by the method `self.force_extension_pair()`.
        xlim : (float, float), optional
            Set xlim of the axis.
        ylim : (float, float), optional
            Set ylim of the axis.
        title : str
            Set the title of the figure.
        draw : bool
            Redraw, after the the data has been updated.
        bps : int
            If lenght of DNA in basepairs is given, a standard stretching
            curve is plotted.
        """
        if self.fe_figure is None:
            self.init_fe_fig(show=False)

        fig = self.fe_figure
        ax = fig.axes[0]

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        # Get stress/release extension/force data pair
        if fe_pair is not None:
            # Unpack force_extension_pair and set the line data to the unpacked
            # values
            time = False
            if len(fe_pair) == 8:
                time = True

            if time:
                est, fst, sti, erl, frl, rli, tst, trl = fe_pair
            else:
                est, fst, sti, erl, frl, rli = fe_pair
            ax.lines[0].set_data(est, fst)
            ax.lines[1].set_data(erl, frl)
            # Autodetermine the title
            if time:
                title = title or "Force Extension, %s, %s, (%.3f, %.3f) s" \
                                 % (sti[0], sti[1], tst[0], trl[-1])
            else:
                title = title or "Force Extension, %s, %s" % (sti[0], sti[1])

        # Calculate force extension of a dna with a known length and plot it
        if bps:
            x, F = dna.force_extension(bps=bps)
            ax.lines[2].set_data(x*1e9, F*1e12)
        else:
            ax.lines[2].set_data([0], [0])

        if title is not None:
            ax.set_title(title)

        if draw:
            fig.canvas.draw()

    def plot_force_extension(self, i, autolimit=True, xlim=None, ylim=None,
                             bps=None, draw=True, **kwargs):
        """
        Plot the force extension data with index `i` (see method
        `self.force_extension_pairs()`) on self.fe_figure.

        Parameters
        ----------
        i : int
            Index of force extension pair. See method
            `self.force_extension_pair()`.
        autolimit : bool, optional
            Automatically detect xlim and/or ylim from min/max values of force
            and extension. The limits are only detected, if `xlim` and/or
            `ylim` are not supplied.
        xlim : (float, float), optional
            Xlimit of force extension axis.
        ylim : (float, float), optional
            Ylimit of force extension axis.
        bps : int, optional
            Base pairs of standard force extension curve of a dsDNA to be
            plotted along with the experimental data.
        draw : bool, optional
            Update the figure, after data has been plotted.
        **kwargs : dict, optional
            Keyword arguments passed to `self.force_extension_pair()`.
        """
        fe_pair = self.force_extension_pair(i, **kwargs)
        e = np.r_[fe_pair[0], fe_pair[3]]
        f = np.r_[fe_pair[1], fe_pair[4]]
        time = kwargs.get('time', False)

        if autolimit:
            xlim, ylim = self.autolimit(e=e, f=f, xlim=xlim, ylim=ylim,
                                        set_limits=False)
        if time:
            title = "Force Extension %.3i, %s, %s, (%.3f, %.3f) s" \
                    % (i, fe_pair[2][0], fe_pair[2][1], fe_pair[6][0],
                       fe_pair[7][-1])
        else:
            title = "Force Extension %.3i, %s, %s" % (i, fe_pair[2][0],
                                                      fe_pair[2][1])

        self.update_force_extension_plot(fe_pair=fe_pair, xlim=xlim, ylim=ylim,
                                         title=title, draw=draw, bps=bps)

    def plot_force_extensions(self, autolimit=True, xlim=None, ylim=None,
                              bps=None, draw=True):
        """
        Plot force/extension stress/release pairs on an axis of a figure.

        Parameters
        ----------
        autolimit : bool, optional
            Automatically detect xlim and/or ylim from min/max values of force
            and extension. The limits are only detected, if `xlim` and/or
            `ylim` are not supplied.
        xlim : (float, float), optional
            Xlimit of force extension axis.
        ylim : (float, float), optional
            Ylimit of force extension axis.
        bps : int, optional
            Base pairs of standard force extension curve of a dsDNA to be
            plotted along with the experimental data.
        draw : bool, optional
            Update the figure, after data has been plotted.

        Yields
        ------
        matplotlib.figure.Figure
            Yield a reference to the figure, after the next force extension
            force/extension stress/release data pair has been plotted.
        """
        # Get xlim and ylim according to min/max values of the force extension
        # pairs to be plotted
        if autolimit:
            xlim, ylim = self.autolimit(xlim=xlim, ylim=ylim, set_limits=False)

        # Get force extension pairs generator
        fe_pairs = self.force_extension_pairs()

        # Plot all stress/release extension/force data pairs
        for idx, fe_pair in enumerate(fe_pairs):
            title = "Force Extension %.3i, %s, %s" % (idx, fe_pair[2][0],
                                                      fe_pair[2][1])
            self.update_force_extension_plot(fe_pair=fe_pair, xlim=xlim,
                                             ylim=ylim, title=title, draw=draw,
                                             bps=bps)
            yield self.fe_figure

    def show_force_extension_plots(self, **kwargs):
        """
        Display a slider to show the force extensions created by
        `self.plot_force_extension()`.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments passed to `self.plot_force_extension()`.

        Returns
        -------
        ipywidgets.IntSlider
        """
        # Get number of all force extension pairs
        stop = len(self.stress_release_pairs(**kwargs)[0])

        def view_image(i):
            try:
                self.plot_force_extension(i, **kwargs)
            except:
                print('No image found!')

        slider = IntSlider(min=0, max=stop - 1, step=1, value=0,
                           description='Image:')

        return interact(view_image, i=slider)

    def save_force_extension_plots(self, directory=None, file_prefix=None,
                                   file_suffix=None, file_extension='.png',
                                   **kwargs):
        """
        Save all plots created by `plot_force_extensions()`.

        directory : str
            The directory the images to be displayed are located in.
        file_prefix : str
            Display only the files beginning with `prefix`.
        file_suffix : str
            Display only the files ending with `suffix`.
        file_extension : str, optional
            The extension of the images that should be displayed. Default is
            '.png'.
        figure : matplotlib figure, optional
            A reference to a figure that should be used to plot the force
            extension pairs. If no figure is given, a new one is automatically
            created.
        **kwargs
            Parameters passed to the method `self.plot_force_extensions()`.
        """
        kwargs.pop('draw', None)
        # Create generator for all force/extension stress/release pairs
        figures = self.plot_force_extensions(draw=False, **kwargs)

        # Save all figures
        evaluate.save_figures(figures, directory=directory,
                              file_prefix=file_prefix, file_suffix=file_suffix,
                              file_extension=file_extension)

        # Redraw the figure, after the last one has been saved
        self.fe_figure.canvas.draw()

    @property
    def displacementXYZ(self):
        """
        Displacement in µm with height dependent calibration factors for X, Y
        and Z.
        """
        data = self.get_data(traces=['psdXYZ', 'positionZ'])
        psdXYZ = data[:, 0:3]
        positionZ = data[:, [3]]
        displacementXYZ = self.calibration.displacement(psdXYZ,
                                                        positionZ=positionZ)
        return displacementXYZ

    @property
    def forceXYZ(self):
        """
        Force in nN, that is acting on the tether
        """
        positionZ = self.get_data(traces='positionZ')
        forceXY_Z = self.calibration.force(self.displacementXYZ,
                                           positionZ=positionZ)
        return forceXYZ(forceXY_Z)

    def _force(self, samples=None, twoD=False):
        """
        Magnitude of the force in nN acting on the tethered molecule (1D
        numpy.ndarray).
        """
        # Get force (in a fast way)
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)

        psdXYZ = data[:, 0:3]
        positionXY = data[:, 3:5]
        positionZ = data[:, [5]]

        displacementXYZ \
            = self.calibration.displacement(psdXYZ, positionZ=positionZ)
        # 2D or 3D calculation of the distance in Z
        if twoD:
            displacementXYZ[:, Z] = 0.0
        forceXY_Z = self.calibration.force(displacementXYZ,
                                           positionZ=positionZ)

        fXYZ = forceXYZ(forceXY_Z)
        f = force(fXYZ, positionXY)

        return f

    @property
    def force(self):
        """
        Magnitude of the force in nN determined in 3D, acting on the tethered
        molecule.
        """
        return self._force()

    @property
    def distanceXYZ(self):
        """
        Distance of the attachment point to the bead center for all 3 axes.
        """
        # µm, point of attachment of DNA
        positionXYZ = self.get_data(traces='positionXYZ')
        return distanceXYZ(positionXYZ, self.displacementXYZ,
                           self.calibration.radius,
                           self.calibration.focalshift)

    @property
    def distance(self):
        """
        Distance of the attachment point to the bead center.
        """
        positionXY = self.get_data(traces='positionXY')
        return distance(self.distanceXYZ, positionXY)

    def _extension(self, samples=None, twoD=False):
        """
        Extension in µm of the tethered molecule (1D numpy.ndarray).
        """
        # Get extension (in a fast way)
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)

        psdXYZ = data[:, 0:3]
        positionXYZ = data[:, 3:6]
        positionXY = positionXYZ[:, 0:2]
        positionZ = data[:, [5]]

        displacementXYZ \
            = self.calibration.displacement(psdXYZ, positionZ=positionZ)
        # 2D or 3D calculation of the distance in Z
        if twoD:
            displacementXYZ[:, Z] = 0.0
        distXYZ = distanceXYZ(positionXYZ, displacementXYZ,
                              self.calibration.radius,
                              self.calibration.focalshift)

        dist = distance(distXYZ, positionXY)
        e = extension(dist, self.calibration.radius)

        return e

    @property
    def extension(self):
        """
        Extension of the tethered molecule in µm.
        """
        return self._extension()

    def _force_extension(self, samples=None, twoD=False):
        """
        Extension (µm, first column) of and force (nN, second column) acting
        on the tethered molecule (2D numpy.ndarray).
        """
        # Get extension and force (in a fast way)
        data = self.get_data(traces=['psdXYZ', 'positionXYZ'], samples=samples)

        psdXYZ = data[:, 0:3]
        positionXYZ = data[:, 3:6]
        positionXY = positionXYZ[:, 0:2]
        positionZ = data[:, [5]]

        displacementXYZ \
            = self.calibration.displacement(psdXYZ, positionZ=positionZ)
        # 2D or 3D calculation of the distance in Z
        if twoD:
            displacementXYZ[:, Z] = 0.0
        distXYZ = distanceXYZ(positionXYZ, displacementXYZ,
                              self.calibration.radius,
                              self.calibration.focalshift)

        dist = distance(distXYZ, positionXY)
        forceXY_Z = self.calibration.force(displacementXYZ,
                                           positionZ=positionZ)

        fXYZ = forceXYZ(forceXY_Z)

        e = extension(dist, self.calibration.radius)
        f = force(fXYZ, positionXY)

        return np.c_[e, f]

    @property
    def force_extension(self):
        """
        Extension (µm, first column) of and force (nN, second column) acting
        on the tethered molecule (2D numpy.ndarray).
        """
        return self._force_extension()

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


def forceXYZ(forceXY_Z, copy=True):
    """
    Force in nN, that is acting on the tether
    """
    if copy:
        forceXYZ = forceXY_Z.copy()
    else:
        forceXYZ = forceXY_Z
    # stressing on the bead (negative displacement) corresponds to a positive
    # force in Z
    forceXYZ[:, Z] = - 1.0 * forceXY_Z[:, Z]
    return forceXYZ


def force(forceXYZ, positionXY):
    """
    Returns force.

    Parameters
    ----------
    forceXYZ : 2D numpy.ndarray of type float
        forceXYZ.shape[1] can consist of either 3 (XYZ) or 2 (XY) axes
    """
    # sign of forceXY depends on positionXY, important for noise of force
    # around +/- 0
    # forceZ negative/positive, irrespective of positionZ!
    signF = np.sign(forceXYZ)
    signF[:, XY] = np.sign(forceXYZ[:, XY]) * np.sign(positionXY)

    # square the forces and account for the signs
    force_sq = forceXYZ**2 * signF
    forceSUM = np.sum(force_sq, axis=1)
    force = np.sqrt(np.abs(forceSUM)) * np.sign(forceSUM)
    return force


def distanceXYZ(positionXYZ, displacementXYZ, radius=0.0, focalshift=1.0,
                clip_Z=True):
    """
    Distance of the attachment point to the bead center for all 3 axes.
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
    # distance, point of attachment of DNA
    # displacement, displacement of bead out of trap center
    # radius
    distanceXYZ = positionXYZ.copy()
    # distance from attachment point to center of bead
    # attachmentXY - displacementXY
    distanceXYZ[:, 0:2] -= displacementXYZ[:, 0:2]

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
                                + displacementXYZ[idx_free, 2])
    distanceXYZ[idx_touch, 2] = (- positionXYZ[idx_touch, 2]
                                 + radius
                                 # distanceZ + radius + displacementZ
                                 + displacementXYZ[idx_touch, 2])

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


def distance(distanceXYZ, positionXY):
    """
    Calculate the distance of the attachment point to the bead center.

    Parameters
    ----------
    distanceXYZ : 2D numpy.ndarray of type float
        distanceXYZ.shape[1] can consist of either 3 (XYZ) or 2 (XY) axes
    """
    # sign of distanceXY depends on positionXY, important for noise of dist
    # around +/- 0
    # distanceZ negative/positive, irrespective of positionZ!
    signD = np.sign(distanceXYZ)
    signD[:, XY] = np.sign(distanceXYZ[:, XY]) * np.sign(positionXY)

    # square the distances and account for the signs
    distance_sq = distanceXYZ**2 * signD
    distSUM = np.sum(distance_sq, axis=1)
    return np.sqrt(np.abs(distSUM)) * np.sign(distSUM)


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


def _force_extension_figure(grid=True):
    """
    Create a figure to plot force/extension stress/release pair(s).

    Parameters
    ----------
    grid : bool
        Switch on/of the grid of the axis.

    Returns
    -------
    matplotlib.figure.Figure
        Return a reference to the figure.
    """
    # Create figure and axis for plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Set labels
    ax.set_xlabel("Extension (nm)")
    ax.set_ylabel("Force (pN)")

    # Switch on/off grid
    ax.grid(grid)

    # Create lines for stress/release extension/force data pairs
    # line_stress =
    ax.plot([0], [0], 'm.', ms=1.0)
    # line_release =
    ax.plot([0], [0], 'c.', ms=1.0)
    # standard curve for DNA
    ax.plot([0], [0], 'r')

    # Plot stress/release extension/force data pairs
    title = "Force Extension"
    ax.set_title(title)

    return fig
