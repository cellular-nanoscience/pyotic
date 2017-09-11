# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:58:54 2016

@author: Tobias Jachowski
"""
import numpy as np

from pyoti.modification.modification import Modification
from pyoti import helpers as hp
from pyoti import traces as tc
from pyoti.picklable import unboundfunction


class GenericMod(Modification):
    """
    Generic modification that modifies a set of chosen traces.
    """
    def __init__(self, modify=None, recalculate=None, mod_params=None,
                 print_info=None, **kwargs):
        # modify: function, that modifies the data
        # recalculate: function, which (in automatic mode) is called before
        #               _modify(). Can be used to automatically calculate
        #               mod_params, which are calculated with data obtained
        #               from view_based, i.e. data obtained by
        #               _get_data_based(). (Only called, if self.updated
        #               is False. Self.updated is set to False, upon any change
        #               of the view_based). The mod_params can be used by
        #               _modify().
        # kwargs:
        # traces_apply=None,
        # view_apply=None, # obligatory
        # view_based=None,
        # automatic_switch=False, Option to disable automatic call of
        #           _recalculate() before the call of _modify()
        # datapoints=-1, if>0 creates an iattribute, which is used for
        #                     decimate/average calculation in functions
        #                   _get_data(), _get_data_apply, and _get_data_based()
        #                  and binnumber determination in calculate_bin_means()

        if modify is None:
            raise TypeError("GenericMod missing required positional argument "
                            "`modify`.")

        # initialize _mod_params before call of super().__init__
        # to make sure __getattr__() doesn't end in a infinite
        # recursion loop
        self._mod_params = []

        super().__init__(**kwargs)

        # register widgets for modifications of traces
        for trace in self.traces_apply:
            self._mod_params.append(trace)
            key = self._key(trace)
            description = ''.join(('Parameter for trace ', tc.label(trace)))
            self.add_iattribute(key, description=description, value=0.0)

        if mod_params is not None:
            mod_params = hp.listify(mod_params)
            for name in mod_params:
                self._mod_params.append(name)
                key = self._key(name)
                description = ''.join(('Parameter ', name))
                self.add_iattribute(key, description=description, value=0.0)

        self._mod = unboundfunction(modify)

        if recalculate is not None:
            self._recalc = unboundfunction(recalculate)

        if print_info is not None:
            self._print_inf = unboundfunction(print_info)

    def _modify(self, data, samples, data_traces, data_index, mod_index):
        return self._mod(self, data, samples, data_traces, data_index,
                         mod_index)

    def _recalculate(self):
        if hasattr(self, '_recalc'):
            self._recalc(self)

    def _print_info(self):
        if hasattr(self, '_print_inf'):
            self._print_inf()

    def _key(self, trace):
        return ''.join(('mod_param_', trace))

    def get_mod_params(self, names=None):
        if names is None:
            # trace should have the same dimension as traces_apply
            names = self._mod_params
        else:
            names = hp.listify(names)

        return np.array([self.iattributes[self._key(name)] for name in names])

    def set_mod_params(self, values, names=None, leave_automatic=False):
        if names is None:
            # trace should have the same dimension as mod_params
            names = self._mod_params
        else:
            names = hp.listify(names)
        values = hp.listify(values)

        for name, value in zip(names, values):
            key = self._key(name)
            self.iattributes.set_value(key, value,
                                       leave_automatic=leave_automatic)

    def __getattr__(self, name):
        """
        Allow attributes to be used as params selections for set/get_params
        """
        if name in self._mod_params:
            key = self._key(name)
        # if key in self.iattributes:
            return self.iattributes[key]
        else:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        # key = self._key(name)
        if not self._p_setattr(name, value):
            if '_mod_params' in self.__dict__ and name in self._mod_params:
                self.set_mod_params(value, name)
            else:
                super().__setattr__(name, value)

    mod_params = property(get_mod_params, set_mod_params)
