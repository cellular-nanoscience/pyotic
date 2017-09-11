# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:33:22 2016

@author: Tobias Jachowski
"""
import matplotlib.pyplot as plt
import numpy as np

from pyoti.modification.modification import Modification, GraphicalMod
from pyoti import traces as tc
from pyoti.evaluate import signal as sn


class IAttachment(GraphicalMod):
    """
    Subclass of Attachment that provides graphical interfaces to adjust the fit
    parameters.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Define some transient parameters needed for the graphical fitting of
        # the attachment plateaus [excited_psd, excited_position] for left and
        # right stress views
        self.rightposition = None
        self.leftposition = None
        self.rightpsd = None
        self.leftpsd = None

        # Dimensions of the buttons to adjust the plateaus:
        self.left = 0.13
        self.bottom = 0.79
        self.width = 0.0625
        self.height = 0.046875
        self.bspace = 0.01
        self.lspace = 0.01333

        self._lines = {}
        self._ax = None
        self.supertitle = None
        self._button = {}

        # create some properties of actions/corresponding buttons
        action  = [ 'left',  'right', 'lleft', 'rright', 'up',    'down'  ]
        label   = [ '<',     '>',     '<<',    '>>',     'up',    'down'  ]
        offsetP = [  0.0,     0.0,     0.0,     0.0,     -0.0025,  0.0025 ]
        offsetS = [  0.0025, -0.0025,  0.025,  -0.025,    0.0,     0.0    ]
        row     = [  0,       0,       1,       1,        1,       0      ]
        column  = [  0,       2,       0,       2,        1,       1      ]
        self.action  = action
        self.label   = dict(list(zip(action, label)))
        self.offsetP = dict(list(zip(action, offsetP)))
        self.offsetS = dict(list(zip(action, offsetS)))
        self.row     = dict(list(zip(action, row)))
        self.column  = dict(list(zip(action, column)))

    def _figure(self):
        """
        Initialize and show an interactive plot for adjusting the attachment.
        Adjust the plateau correction parameters interactively and set the
        modification accordingly.
        The plot is stored in self.figure.
        """
        # create new figure and axes for adjusting the plateaus
        figure, ax = plt.subplots(1, 1, sharex=True, sharey=True)
        self._ax = ax

        # create buttons for interactive correction of plateaus and assign
        # correct functions
        # see http://math.andrej.com/2009/04/09/pythons-lambda-is-broken/ for
        # explanation of weird function assignment
        ax_button = {}
        for ac in self.action:
            ax_button[ac] = figure.add_axes([self.column[ac] *
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

        # create lines to plot the plateaus
        self._lines['left'] = ax.plot([0], [0], 'r', alpha=0.75)[0]
        self._lines['right'] = ax.plot([0], [0], 'g', alpha=0.75)[0]
        ax.ticklabel_format(useOffset=False)
        ax.grid(True)

        self.supertitle = figure.suptitle("Adjust plateaus to make them "
                                          "overlap")

        return figure

    def _update_fig(self, **kwargs):
        """
        Update the plot
        """
        # recalculate the plateaus with the new offset and scaling values
        self.calculate_plateaus()

        # plot data
        self._lines['left'].set_data(self.leftposition, self.leftpsd)
        self._lines['right'].set_data(self.rightposition, self.rightpsd)
        excited_psd = self.modification.traces_apply[0]
        excited_position = self.modification.traces_apply[1]
        self._ax.set_xlabel(tc.label(excited_position))
        self._ax.set_ylabel(tc.label(excited_psd))

        # recompute ax.dataLim
        self._ax.relim()
        # update ax.viewLim using new dataLim
        self._ax.autoscale_view()

    def _pre_close_fig(self):
        # Store attachment fit plot for ducumentation
        self.supertitle.set_text('Adjusted  plateaus')
        self._lines.clear()
        self._ax = None
        self._button.clear()
        self.supertitle = None

    def _adjust_plateaus(self, action):
        """
        Adjusts the attachment (offset of excited_position) and the the scaling
        factor to correct for differences of left and right DNA overstretching
        plateaus.
        It is interactively called from the data plot (see below) and updates
        the plot accordingly.
        """
        # change offset and scaling for plateaus
        self.modification.iattributes.offsetStage += self.offsetS[action]
        self.modification.iattributes.offsetPsd += self.offsetP[action]

        self.update_fig()

    def calculate_plateaus(self):
        """
        Calculate the plateaus according to the offsets and the scaling of
        data_based.
        """
        # determine the excited axis of position and psd signal
        ex = self.modification._excited()
        excited_psd = self.modification._NAME['psd'][ex]
        excited_position = self.modification._NAME['position'][ex]
        self.modification.traces_apply = [excited_psd, excited_position]

        # recalculate data for plotting
        # data_based: [excited_psd, excited_position]
        data_based = self.modification._get_data_based(
            traces=self.modification.traces_apply, window=False, decimate=True)

        # subtract offsets
        data_based[:, 0] -= self.modification.iattributes.offsetPsd
        data_based[:, 1] -= self.modification.iattributes.offsetStage

        # determine left/right stress regions of DNA
        signal = data_based[:, 1]  # [excited_position]
        resolution = self.modification.view_based.samplingrate \
            / self.modification.decimate
        minima, maxima = sn.get_extrema(signal, resolution)
        rightstress, _, leftstress = \
            sn.get_sections(signal, minima, maxima)[1][0:3]

        # invert the signals of left stress cycle
        data_based[leftstress] *= -1

        # set plateau data arrays
        self.rightposition = data_based[rightstress, 1]
        self.leftposition = data_based[leftstress, 1]
        self.rightpsd = data_based[rightstress, 0]
        self.leftpsd = data_based[leftstress, 0]


class Attachment(Modification):
    """
    Determine attachment point of DNA and scaling of lateral PSD
    """
    GRAPHICALMOD = IAttachment

    def __init__(self, db_update=False, **kwargs):
        super().__init__(datapoints=25000, **kwargs)

        # determine the excited axis of position and psd signal
        if not db_update:
            ex = self._excited()
            excited_psd = self._NAME['psd'][ex]
            excited_position = self._NAME['position'][ex]
            self.traces_apply = [excited_psd, excited_position]

        # Define parameters that are used to calculate the modification
        # offset of PSD relative to trap center position of bead
        self.add_iattribute('offsetPsd', description='Offset PSD (V)',
                            value=0.0)
        # offset of the position relative to the attachment point of the DNA
        self.add_iattribute('offsetStage', description='Offset position (V)',
                            value=0.0)

    def _print_info(self):
        print("    Excited axis is: %s" % self.traces_apply[1])

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        # correct attachment point of DNA: adjust excited_position
        # (set it to 0 where DNA is attached)
        # correct for offset of excited_psd
        offset = np.array([self.iattributes.offsetPsd,  # excited_psd,
                           self.iattributes.offsetStage])  # excited_position

        data[:, data_index] -= offset[np.newaxis, mod_index]

        return data


# The following is only to update to database version 0.8.0
class GAttachment(Attachment):
    pass
