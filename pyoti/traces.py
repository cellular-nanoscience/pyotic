# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:39:13 2015

@author: Tobias Jachowski
"""
import collections
import sys

from . import config as cf
from . import helpers as hp


class Traces(object):

    def __init__(self, cfgfile, **kwargs):
        self._cfgfile = cfgfile
        self.reload()

    def __getattr__(self, name):
        return self.normalize(name)

    def __contains__(self, key):
        return key in self._TRACE_ALIASES

    def reload(self, cfgfile=None):
        cfgfile = cfgfile or self.cfgfile
        cfg = cf.read_cfg_file(cfgfile)
        self._cfgfile = cfgfile

        self._TRACE_ALIASES = {a: [t.strip() for t in traces.split(',')]
                               for a, traces
                               in cfg['aliases'].items()}

        self._COLOR = {t: c
                       for t, c
                       in cfg['color'].items()}

        self._LABEL = {t: l
                       for t, l
                       in cfg['label'].items()}

    def normalize(self, traces):
        """
        Normalize traces returns a list of trace names.
        It looks up alias names for traces and replaces them with the list of
        names. If alias is not found, the input parameters are converted /
        added to the return list.
        """
        traces = hp.listify(traces)

        if isinstance(traces, collections.Iterable):
            _traces = []
            for trace in traces:
                if isinstance(trace, str) and trace in self._TRACE_ALIASES:
                    # lookup for traces as alias, a shorthand notation, or a
                    # combination
                    trace = self._TRACE_ALIASES[trace]
                    _traces.extend(trace)
                else:
                    _traces.append(trace)

            # Remove duplicate entries
            seen = set()
            seen_add = seen.add
            traces = [trace
                      for trace in _traces
                      if not (trace in seen or seen_add(trace))]

        return traces

    def label(self, trace):
        if trace in self._LABEL:
            return self._LABEL[trace]
        elif 'DEFAULT' in self._LABEL:
            return self._LABEL['DEFAULT']
        else:
            return 'Data (arbitrary)'

    def color(self, trace):
        if trace in self._COLOR:
            return self._COLOR[trace]
        elif 'DEFAULT' in self._COLOR:
            return self._COLOR['DEFAULT']
        else:
            return 'black'

    @property
    def cfgfile(self):
        return self._cfgfile

    @cfgfile.setter
    def cfgfile(self, cfgfile):
        self.reload(cfgfile=cfgfile)

    # ### module like functionality ###
    @property
    def __file__(self):
        return __file__

    @property
    def __name__(self):
        return __name__

    @property
    def __package__(self):
        return __package__

cfgfile='traces.cfg'
sys.modules[__name__] = Traces(cfgfile)
