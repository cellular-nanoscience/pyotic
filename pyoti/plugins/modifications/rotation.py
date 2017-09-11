# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:36:44 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np

from pyoti.modification.modification import Modification, GraphicalMod
from pyoti import helpers as hp
from pyoti.evaluate import signal as sn
from pyoti.evaluate import tether as tr


class IRotation(GraphicalMod):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Define some attributes needed for the graphical visualization of the
        # spline

        # plateaus [excited_psd, excited_position] for left and right stress
        # regions
        self.rightextension = None
        self.leftextension = None
        self.rightforce = {}
        self.leftforce = {}

        # Dimensions of the buttons to adjust the plateaus:
        self.left = 0.1
        self.bottom = 0.79
        self.width = 0.0625
        self.height = 0.046875
        self.bspace = 0.01
        self.lspace = 0.01333

        # create some properties of actions/corresponding buttons
        action  = ['upX',   'downX', 'upY',   'downY', 'upZ',   'downZ', 'rotlX', 'rotrX', 'rotlY', 'rotrY', 'rotlZ', 'rotrZ']
        label   = ['upX',   'downX', 'upY',   'downY', 'upZ',   'downZ', 'rotlX', 'rotrX', 'rotlY', 'rotrY', 'rotlZ', 'rotrZ']
        offsetX = [-0.0025,  0.0025,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0   ]
        offsetY = [ 0.0,     0.0,    -0.0025,  0.0025,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0   ]
        offsetZ = [ 0.0,     0.0,     0.0,     0.0,    -0.0025,  0.0025,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0   ]
        rotateX = [ 0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.5,    -0.5,     0.0,     0.0,     0.0,     0.0   ]
        rotateY = [ 0.0,     0.0,     0.0,     0.0,     0.0,    -0.0,     0.0,     0.0,     0.5,    -0.5,     0.0,     0.0   ]
        rotateZ = [ 0.0,     0.0,     0.0,     0.0,     0.0,    -0.0,     0.0,     0.0,     0.0,     0.0,     0.5,    -0.5   ]
        row     = [ 1,       0,      -6,      -7,      -6,      -7,       1,       0,       1,       0,       1,       0     ]
        column  = [ 6,       6,       0,       0,       6,       6,       0,       0,       1,       1,       2,       2     ]
        self.action  = action
        self.label   = dict(list(zip(action, label)))
        self.offsetX = dict(list(zip(action, offsetX)))
        self.offsetY = dict(list(zip(action, offsetY)))
        self.offsetZ = dict(list(zip(action, offsetZ)))
        self.rotateX = dict(list(zip(action, rotateX)))
        self.rotateY = dict(list(zip(action, rotateY)))
        self.rotateZ = dict(list(zip(action, rotateZ)))
        self.row     = dict(list(zip(action, row)))
        self.column  = dict(list(zip(action, column)))

        self._axes = ['3D', 'X', 'Y', 'Z']

        self._lines = {}
        self.supertitle = None
        self._button = {}

    def _figure(self):
        """
        Interactive determination of rotation angle.
        """
        # Adjust the angles for rotation of of the QPD signal

        # create new figure and axes for adjusting the angles
        plot_params = {'wspace': 0.09, 'hspace': 0.08}
        self._set_plot_params(plot_params)
        figure, ax_ = plt.subplots(2, 2, sharex=True, sharey=False)

        ax = dict(list(zip(self._axes, ax_.flatten())))

        # Create buttons for interactive correction of plateaus and assign
        # correct functions
        # See http://math.andrej.com/2009/04/09/pythons-lambda-is-broken/ for
        # explanation of weird function assignment
        ax_button = {}
        for ac in self.action:
            ax_button[ac] \
                = figure.add_axes([self.column[ac] *
                                   (self.lspace + self.width) +
                                   self.left, self.row[ac] *
                                   (self.bspace + self.height) +
                                   self.bottom,
                                   self.width,
                                   self.height])
            self._button[ac] = plt.Button(ax_button[ac], self.label[ac])

            def ap(event, ac=ac):
                self._adjust_plateaus(ac)
            # connect button to action, accordingly
            self._button[ac].on_clicked(ap)

        # create lines to plot the data
        for plot in self._axes:
            self._lines['left' + plot] = ax[plot].plot([0], [0], 'r',
                                                       alpha=0.75)[0]
            self._lines['right' + plot] = ax[plot].plot([0], [0], 'g',
                                                        alpha=0.75)[0]
            ax[plot].ticklabel_format(useOffset=False)
            ax[plot].grid(True)
        ax['3D'].set_ylabel('Force (pN)')
        ax['Y'].set_ylabel('Force (pN)')
        ax['Y'].set_xlabel('Extension (nm)')
        ax['Z'].set_xlabel('Extension (nm)')

        self.supertitle = figure.suptitle("Adjust plateaus to make them "
                                          "overlap")

        return figure

    def _update_fig(self, **kwargs):
        """
        Update the plot
        """
        # prepare data for interactive plot
        self.calculate_plateaus()

        for plot in self._axes:
            l_line = self._lines['left' + plot]
            r_line = self._lines['right' + plot]

            l_line.set_data(self.leftextension * 1000,
                            self.leftforce[plot] * 1000)
            r_line.set_data(self.rightextension * 1000,
                            self.rightforce[plot] * 1000)
            # recompute ax.dataLim
            l_line.axes.relim()
            # update ax.viewLim using new dataLim
            l_line.axes.autoscale_view()

    def _pre_close_fig(self):
        # Store attachment fit plot for ducumentation
        self.supertitle.set_text('Adjusted  plateaus')
        self._lines.clear()
        self._button.clear()
        self.supertitle = None

    def _adjust_plateaus(self, action):
        """
        Adjusts the attachment (offset of excited_position) and the the scaling
        factor to correct for differences of left and right DNA overstretching
        plateaus. It is interactively called from the data plot (see below) and
        updates the plot accordingly.
        """
        # change offset and scaling for plateaus
        self.modification.iattributes.offsetPsdX += self.offsetX[action]
        self.modification.iattributes.offsetPsdY += self.offsetY[action]
        self.modification.iattributes.offsetPsdZ += self.offsetZ[action]

        self.modification.iattributes.angleX += self.rotateX[action]
        self.modification.iattributes.angleY += self.rotateY[action]
        self.modification.iattributes.angleZ += self.rotateZ[action]

        # recalculate the plateaus with the new offset and scaling values
        self.update_fig()

    def calculate_plateaus(self):
        """
        Calculate the plateaus according to the offsets and the angles of
        data_based.
        """
        # recalculate data for plotting
        traces = ['psdX', 'psdY', 'psdZ', 'positionX', 'positionY',
                  'positionZ']
        data_based = self.modification._get_data_based(
            traces=traces, window=False, decimate=True)

        # subtract offsets
        data_based[:, 0] -= self.modification.iattributes.offsetPsdX
        data_based[:, 1] -= self.modification.iattributes.offsetPsdY
        data_based[:, 2] -= self.modification.iattributes.offsetPsdZ

        # set positionZ and calibration, needed by
        # self.rotate -> self.rot_factor
        # make sure that both position signal and calibration are taken form
        # self.view_apply
        positionZ = data_based[:, hp.slicify(5)]
        calibration = self.modification.view_based.calibration

        # rotate psd
        data_based[:, 0], data_based[:, 1], data_based[:, 2] \
            = self.modification.rotate(
                data_based[:, 0], data_based[:, 1],
                data_based[:, 2], positionZ, calibration)

        # set some variable names for easy access of data_based
        psdXYZ = data_based[:, hp.slicify([0, 1, 2])]
        positionXYZ = data_based[:, hp.slicify([3, 4, 5])]
        positionXY = data_based[:, hp.slicify([3, 4])]
        positionZ = data_based[:, hp.slicify(5)]

        # calculate extension
        displacementXYZ = calibration.displacement(psdXYZ, positionZ=positionZ)
        distanceXYZ = tr.distanceXYZ(positionXYZ, displacementXYZ,
                                     calibration.radius,
                                     calibration.focalshift)
        distance = tr.distance(distanceXYZ, positionXY)
        extension = tr.extension(distance, calibration.radius)

        # calculate force
        force3D = calibration.force(displacementXYZ, positionZ=positionZ)
        force = tr.force(tr.forceXYZ(force3D), positionXY)
        force = {'3D': force,
                 'X': force3D[:, 0],
                 'Y': force3D[:, 1],
                 'Z': force3D[:, 2]}

        # determine regions where DNA is stretched to the right and left side
        ex = self.modification._excited()
        excited_position = self.modification._NAME['position'][ex]
        positionEl = traces.index(excited_position)
        signal = data_based[:, positionEl]  # [positionE]
        resolution = self.modification.view_based.samplingrate \
            / self.modification.decimate
        minima, maxima = sn.get_extrema(signal, resolution)
        rightstress, _, leftstress \
            = sn.get_sections(signal, minima, maxima)[1][0:3]

        # set plateau data arrays
        self.rightextension = extension[rightstress]
        self.leftextension = extension[leftstress]
        for plot in self._axes:
            self.rightforce[plot] = force[plot][rightstress]
            self.leftforce[plot] = force[plot][leftstress]


class Rotation(Modification):
    GRAPHICALMOD = IRotation

    def __init__(self, **kwargs):
        traces_apply = ['psdX', 'psdY', 'psdZ']
        super().__init__(datapoints=12000, traces_apply=traces_apply, **kwargs)

        # Define parameters that are used to calculate the modification
        # the angles the ellipsoid has to be rotated
        # rotation around X
        self.add_iattribute('angleX', description='Angle about X (deg)',
                            value=0.0)
        self.add_iattribute('angleY', description='Angle about Y (deg)',
                            value=0.0)
        self.add_iattribute('angleZ', description='Angle about Z (deg)',
                            value=0.0)
        # offset of PSD relative to trap center position of bead
        self.add_iattribute('offsetPsdX', description='Offset PSD X (V)',
                            value=0.0)
        self.add_iattribute('offsetPsdY', description='Offset PSD Y (V)',
                            value=0.0)
        self.add_iattribute('offsetPsdZ', description='Offset PSD Z (V)',
                            value=0.0)

        # Parameters for rotation matrix calculation
        self.rotation_method = 'nm'  # 'pN', 'nm' or 'V'

    def _print_info(self):
        print(("    Rotation is in '%s' space" % self.rotation_method))

    def _modify(self, data, samples, data_traces, data_index, mod_index):

        # Get data that is needed for the rotation modification, but is not
        # contained in the data array, that is requested to be getting modified
        # (`data_traces` in `data`)
        needed_traces = self.traces_apply  # traces_apply
        needed_traces.append('positionZ')  # psdX, psdY, psdZ, positionZ
        # calculate missing traces
        extra_traces = hp.missing_elements(needed_traces, data_traces)
        extra_data = self._get_data_apply(samples=samples, traces=extra_traces,
                                          copy=False)

        # combine assigned data and traces with extra fetched data and traces
        traces_tuple = (data_traces, extra_traces)
        data_tuple = (data, extra_data)

        # function to easily get data for a trace from different data
        def get_target_data(target_trace, traces_tuple, data_tuple):
            target_data = None
            for i, traces in enumerate(traces_tuple):
                if target_trace in traces:
                    index = traces.index(target_trace)
                    target_data = data_tuple[i][:, index]
            return target_data

        # get psdX, psdY and psdZ
        l_data = {}
        for trace in self.traces_apply:
            l_data[trace] = get_target_data(trace, traces_tuple, data_tuple)

        # correct for offset of psds
        x = self.iattributes.offsetPsdX
        y = self.iattributes.offsetPsdY
        z = self.iattributes.offsetPsdZ
        if x != 0.0:
            l_data['psdX'] = l_data['psdX'] - x
        if y != 0.0:
            l_data['psdY'] = l_data['psdY'] - y
        if z != 0.0:
            l_data['psdZ'] = l_data['psdZ'] - z

        # positionZ and calibration, needed by
        # self.rotate -> self.rot_factor
        # make sure that both position signal and calibration are taken form
        # self.view_apply
        positionZ = get_target_data('positionZ', traces_tuple, data_tuple)
        calibration = self.view_apply.calibration

        # rotate the data
        l_data['psdX'], l_data['psdY'], l_data['psdZ'] \
            = self.rotate(l_data['psdX'], l_data['psdY'], l_data['psdZ'],
                          positionZ, calibration)

        for trace in self.traces_apply:
            if trace in data_traces:
                index = data_traces.index(trace)
                data[:, index] = l_data[trace]

        return data

    def rotate(self, data_x, data_y, data_z, positionZ, calibration):
        ax = self.iattributes.angleX
        ay = self.iattributes.angleY
        az = self.iattributes.angleZ

        if ax == 0.0 and ay == 0.0 and az == 0.0:
            return (data_x, data_y, data_z)

        x = ax * np.pi / 180.0
        y = ay * np.pi / 180.0
        z = az * np.pi / 180.0

        cf = self.calibration_factor(positionZ, calibration)
        rf = self.rot_factor

        # https://www.siggraph.org/education/materials/HyperGraph/modeling/mod_tran/3drota.htm
        # angle about axis ...
        x_o = data_x
        y_o = data_y
        z_o = data_z

        # right handed coordinate system
        # rotation is counterclockwise about the axis coming out of the image
        # plane
        # z should always be <= 0 (negative), because the bead is stressed down
        #
        # rotate about x
        #
        #    z |
        #      |
        #      |
        #      |_______ y
        #      /
        #     /
        #  x /
        #
        if x == 0.0:
            y_2 = y_o
            z_2 = z_o
        else:
            cos_x = np.cos(x)
            sin_x = np.sin(x)
            y_2 = y_o * cos_x - z_o * rf(2, 1, cf) * sin_x
            z_2 = z_o * cos_x + y_o * rf(1, 2, cf) * sin_x

        # rotate about y
        #
        #    x |
        #      |
        #      |
        #      |_______ z
        #      /
        #     /
        #  y /
        #
        if y == 0.0:
            z_o = z_2
            x_2 = x_o
        else:
            cos_y = np.cos(y)
            sin_y = np.sin(y)
            z_o = z_2 * cos_y - x_o * rf(0, 2, cf) * sin_y
            x_2 = x_o * cos_y + z_2 * rf(2, 0, cf) * sin_y

        # rotate about z
        #
        #    y |
        #      |
        #      |
        #      |_______ x
        #      /
        #     /
        #  z /
        #
        if z == 0.0:
            x_o = x_2
            y_o = y_2
        else:
            cos_z = np.cos(z)
            sin_z = np.sin(z)
            x_o = x_2 * cos_z - y_2 * rf(1, 0, cf) * sin_z
            y_o = y_2 * cos_z + x_2 * rf(0, 1, cf) * sin_z

        return (x_o, y_o, z_o)

    def calibration_factor(self, positionZ, calibration):
        # position signal, to calculate height dependent calibration factors
        if positionZ.ndim < 2:
            positionZ = positionZ[:, np.newaxis]

        # rotation in nm space
        if self.rotation_method == 'nm':
            # beta in nm/mV, displacement sensitivity
            beta = calibration.beta(positionZ)[:, 0:3]
            factor = beta
        # rotation in pN space
        elif self.rotation_method == 'pN':
            # beta in nm/mV, displacement sensitivity
            beta = calibration.beta(positionZ)[:, 0:3]
            # pN/nm, stiffness
            kappa = calibration.kappa(positionZ)[:, 0:3]
            factor = kappa * beta
        else:
            # self.rotation_method == 'V'  # rotation in V space
            factor = 1.0

        return factor

    def rot_factor(self, a, b, calibration_factor):
        """
        Calculates a weighting factor for rotating the 3D QPD signal in V via
        the signal in pN (stiffness * displacement_sensitivity), nm
        (displacement_sensitivity), or V space.
        (If you rotate a vector in pN space, the absolute value (force) should
        stay the same, regardless of the axis. A rotation in V space would lead
        to bias, due to the different displacement sensitivities and
        stiffnesses of all three axes.)
        To calculate the factors, this method uses the calibration of the
        view, this modification is applied to.
        """
        cf = calibration_factor
        if self.rotation_method == 'nm' or self.rotation_method == 'pN':
            factor = cf[:, a] / cf[:, b]
        else:
            factor = 1.0

        return factor

    def Rxyz(self, positionZ, calibration, x=0.0, y=0.0, z=0.0):
        """
        Create rotationmatrix along x, y and z in degrees.
        Rotation for data with N samples and 3 dimensions (XYZ):
            R = Rxyz(x,y,z)
            data_rot = np.dot(data, R.T)

        Rxyz takes the mean() of the rotation factors, calculated by
        self.rot_factor()!
        """

        # angle about axis ...
        x = x * np.pi / 180.0
        y = y * np.pi / 180.0
        z = z * np.pi / 180.0

        f = self.rot_factor  # calculate weighting factors

        fXY = f(0, 1, positionZ, calibration).mean()
        fXZ = f(0, 2, positionZ, calibration).mean()
        fYX = f(1, 0, positionZ, calibration).mean()
        fYZ = f(1, 2, positionZ, calibration).mean()
        fZX = f(2, 0, positionZ, calibration).mean()
        fZY = f(2, 1, positionZ, calibration).mean()

        Rx = np.matrix([[     1.0,            0.0,            0.0       ],
                        [     0.0,            np.cos(x), -fZY*np.sin(x) ],
                        [     0.0,        fYZ*np.sin(x),      np.cos(x) ]])
        Ry = np.matrix([[     np.cos(y),      0.0,       -fZX*np.sin(y) ],
                        [     0.0,            1.0,            0.0       ],
                        [ fXZ*np.sin(y),      0.0,            np.cos(y) ]])
        Rz = np.matrix([[     np.cos(z), -fYX*np.sin(z),      0.0       ],
                        [ fXY*np.sin(z),      np.cos(z),      0.0       ],
                        [     0.0,            0.0,            1.0       ]])

        return Rx*Ry*Rz


# The following is only to update to database version 0.8.0
class GRotation(Rotation):
    pass
