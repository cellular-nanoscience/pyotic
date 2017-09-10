# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:11:44 2016

@author: Tobias Jachowski
"""
import numpy as np

from pyoti.modification.modification import Modification
from pyoti import traces as tc


class Offset(Modification):
    """
    Modification that corrects the offset of psdX, psdY, and psdZ.
    """
    def __init__(self, **kwargs):
        traces_apply = ['psdX', 'psdY', 'psdZ']
        super().__init__(automatic_switch=True, traces_apply=traces_apply,
                         **kwargs)

        # register widgets for the offsets
        for trace in self.traces_apply:
            key = self._key(trace)
            description = ''.join(('Offset ', tc.label(trace)))
            self.add_iattribute(key, description=description, value=0.0)

    def _recalculate(self):
        # Calculate modification for offset from the selected offset span
        traces = self.traces_apply
        data = self._get_data_based(traces=traces, copy=False)
        # should have the same dimension as traces_apply
        self.set_offset(data.mean(axis=0), leave_automatic=True)

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        data[:, data_index] -= self.offset[np.newaxis, mod_index]
        return data

    def _key(self, trace):
        return ''.join(('offset_', trace))

    def get_offset(self):
        return np.array([self.iattributes[self._key(trace)]
                         for trace in self.traces_apply])

    def set_offset(self, offset, trace=None, leave_automatic=False):
        # trace should have the same dimension as traces_apply
        if trace is None:
            traces = self.traces_apply
        else:
            traces = [trace]

        for trace in traces:
            index = self.lia(trace)
            key = self._key(trace)
            self.iattributes.set_value(key, offset[index],
                                       leave_automatic=leave_automatic)

    offset = property(get_offset, set_offset)
