# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 17:09:09 2016

@author: Tobias Jachowski
"""
import hashlib
import inspect
import numpy as np

from .. import config as cf
from .. import traces as tc
from .region import Region


class Record(Region):
    """
    A Record consists of:
        - a filename, containing the name of the filename
        - a parafile, containing the name of the parafile
        - a calibration object, containing dsurf, focalshift and radius
          (and the height dependent stiffnesses and displacement sensitivities)
        - the data loaded from the filename
        - the parameters loaded from the parafile
        - and some methods for convenient access of the data and para.

    A Record contains some basic modifications, like a static and fixed offset,
    conversion, and inversion factor. The offset, inversion, and conversion
    factors depend on the parameters given by the experimental setup and
    software, only.
    The offset factor should be applied only, if the graphical user interface
    of the OT has a different value than the final value stored in the data.
    The inversion factor should be applied only, if the values e.g. of the
    position are moving in an inverted fashion.
    The conversion factor should be applied only, if the values e.g. of the
    position are stored in Volts. The position values are needed to be in um.

    The modifications are needed to have the same values as used during the
    calibration and need to be set accordingly.

    For future development one should consider to implement the OT control
    software such that whereever possible, Volts are converted to SI units and
    a common convention for the directions (left/right, down/up is neg/pos) is
    followed.
    """
    def __init__(self, datasource, traces, calibration, offset=None,
                 conversion=None, inversion=None, **kwargs):
        """
        Parameters
        ----------
        datasource : DataSource
            The object to retrieve data from. This object has to have the
            function as_array() and the attribue `name`.
        offset: dictionary of trace: value pairs
            value can either be a value or a numpy function (e.g. median)

        """
        super().__init__(max_parents=1, caching=True, **kwargs)

        if datasource is None:
            raise TypeError("Record missing the required positional argument "
                            "'datasource'.")
        self._datasource = datasource
        self._dataident = None

        # Instance of Calibration used to get dsurf, radius and focalshift
        if calibration is None:
            raise TypeError("Record missing the required positional argument "
                            "'calibraion'.")
        self.calibration = calibration

        # Initialize self._shape before trying to set offset, inversion, or
        # conversion, which in turn need to (indirectly) access self._shape.
        # ZODB volatile.
        self._v_data_cached = self._raw_data
        raw_data = self._v_data_cached
        self._shape = raw_data.shape

        # Check if there are as many descriptions as traces in the data
        traces = tc.normalize(traces)
        if len(traces) != self.num_traces:
            raise TypeError('Record argument "traces" has {} elements, but '
                            '"datasource" has {} data columns. Please, '
                            'provide the argument "traces" with a list of '
                            'names for every data column, or edit the '
                            '"cfgfile", accordingly'.format(len(traces),
                                                            self.num_traces))
        self._traces = traces

        # Initialize values of offset (0), inversion (1), and conversion (1)
        # for all available traces
        self._offset = np.zeros(self.num_traces)
        self._conversion = np.ones(self.num_traces)
        self._inversion = np.ones(self.num_traces)

        # Modify offset, conversion and inversion according to optional
        # parameters
        try:
            for trace, value in offset.items():
                if isinstance(value, str):
                    trace_idx = self.traces_to_idx(trace)
                    value = getattr(np, value)(raw_data[:, trace_idx], axis=0)

                self.set_offset(trace, value)
        except:
            pass

        try:
            for trace, value in conversion.items():
                self.set_conversion(trace, value)
        except:
            pass

        try:
            for trace, value in inversion.items():
                self.set_inversion(trace, value)
        except:
            pass

        # reset cache according to modified offset, inversion, and conversion.
        # ZODB volatile.
        self._v_data_cached = ((raw_data - self.offset) * self.inversion
                               * self.conversion)

    def _get_data_uncached(self, samples, traces_idx, copy=True):
        """
        Returns the data of Record, according to self.filename.
        Copy is always True for Record, because every time it reads in the data
        it createas a new numpy.ndarray.
        """
        # read in data
        data = ((self._raw_data[samples, traces_idx] - self.offset[traces_idx])
                * self.inversion[traces_idx]
                * self.conversion[traces_idx])
        # TODO: Implement different samplingrates for different traces. Pandas?
        return data

    @property
    def _raw_data(self):
        """
        Reads and returns uncorrected, uncached data (no offset, no inversion,
        no conversion).
        Additionally, sets self._shape according to _raw_data.shape and creates
        or compares previously generated hashdigest of the data.
        """
        data = self.datasource.as_array()
        if not isinstance(data, np.ndarray):
            raise TypeError("The data you try to load is no numpy array!")
        if data.ndim != 2:
            raise ValueError("The data array you try to load does not have 2 "
                             "dimensions!")
        data = data.copy(order='C')
        ident = hashlib.md5(data).hexdigest()
        if self._dataident is None:
            self._dataident = ident
        elif self._dataident != ident:
            raise ValueError("The data you try to load from '%s' has changed "
                             "since the last time. Please, check the "
                             "datasource of the record '%s'."
                             % (self.datasource.name, self.name))

        return data

    @property
    def offset(self):
        return self._offset.copy()

    @property
    def inversion(self):
        return self._inversion.copy()

    @property
    def conversion(self):
        return self._conversion.copy()

    def set_offset(self, trace, offset):
        """
        Set the offset of a selected trace. The offset should be given in units
        of the raw data, this means independent of inversion and conversion.
        """
        trace = self.traces_to_idx(trace)
        self._offset[trace] = offset
        # Inform descendants of change
        self.set_changed()

    def set_inversion(self, trace, invert=True):
        trace = self.traces_to_idx(trace)
        self._inversion[trace] = - int(invert) or 1
        # Inform descendants of change
        self.set_changed()

    def set_conversion(self, trace, conversion):
        trace = self.traces_to_idx(trace)
        self._conversion[trace] = conversion
        # Inform descendants of change
        self.set_changed()

    @property
    def calibration(self):
        return self.parent

    @calibration.setter
    def calibration(self, calibration):
        if calibration is None:
            raise TypeError("Calibration must not be None.")
        self.set_parent(calibration)

    @property
    def datasource(self):
        return self._datasource

    @property
    def samplingrate(self):
        # TODO: Implement different samplingrates for different traces. Pandas?
        return self.datasource.samplingrate

    @property
    def start(self):
        return 0

    @property
    def stop(self):
        return self._shape[0]

    @property
    def num_traces(self):
        return self._shape[1]

    @property
    def traces(self):
        return self._traces.copy()

    @property
    def caching(self):
        return self._caching

    @caching.setter
    def caching(self, caching):
        super(Record, self.__class__).caching.fset(self, caching)
        if caching is False:
            print("You switched of caching for Record.",
                  "This will seriously lessen the performance of pyoti!",
                  "Do this only, if you have weighty reasons.",
                  "Otherwise, you should revert and reset `caching` to True",
                  sep="\n")


class SpecialRecord(Record):
    """
    This class inherits from the superclass Record and can be used as a
    template to modify the default behaviour of how the data is read in.
    Details see property _raw_data().
    """
    def __init__(self, datasource, traces, calibration, offset=None,
                 conversion=None, inversion=None, **kwargs):

        super().__init__(datasource, traces, calibration, offset=offset,
                         conversion=conversion, inversion=inversion, **kwargs)

        # Get corrected raw_data, which by super().__init__() was initialized
        # and stored in self._v_data_cached
        # raw_data = self._v_data_cached

        # Or get uncorrected raw_data
        # raw_data = self._raw_data()

        # Do something to the data for initialization of SpecialRecord

    @property
    def _raw_data(self):
        data = super()._raw_data
        # do something to data ...
        # e.g. calculate positionXYZ of a combination of signals of one or
        # multiple stage, objective, moveable lens, or mirror signals. Keep in
        # mind that positionZ is increasing for increasing distances and
        # decreasing for decreasing distances of the bead to the surface.
        # PositionXY need to have the same sign as psdXY, i.e. if the bead is
        # pulled to the left/right and displaced to the left/right, positionXY
        # and psdXY both need to have the same signs.
        return data


def create_record(calibration,
                  traces=None,
                  name=None,
                  group=None,
                  rc_class=None,
                  ds_class=None,
                  offset=None,
                  conversion=None,
                  inversion=None,
                  cfgfile=None,
                  **kwargs):
    """
    Parameters
    ----------
    **kwargs
        Is used to get parameters for initialization of `rc_class` and
        initialisation of `ds_class`
    """

    # Set default configfile
    cfgfile = cfgfile or 'record.cfg'

    # Read configfile
    cfg = cf.read_cfg_file(cfgfile)

    print(("Creating record '" + name + "':"))

    rc_class = rc_class or cf.get_cfg_class(cfg, sec='record',
                                            std_mod='.region',
                                            std_cls='Record')
    traces = traces or cf.get_cfg_list(cfg, sec='record', opt='traces')

    if not rc_class:
        print("Could not create Record class defined in config file %s"
              % cfgfile)
        return None

    # get parameters specific for record
    record_pars = {}
    for par in inspect.getargspec(rc_class.__init__)[0]:
        if par in kwargs:
            record_pars[par] = kwargs.pop(par)

    # datasource class and parameters
    ds_class = ds_class or cf.get_cfg_class(cfg, sec='datasource',
                                        std_mod='.plugins.datasources.generic',
                                        std_cls='GenericDataFile')

    if not ds_class:
        print("Could not create DataSource class defined in config file %s"
              % cfgfile)
        return None

    # create datasource
    datasource = ds_class(**kwargs)

    # read and convert offset, conversion, and inversion to dict of floats
    # (booleans)
    offset = offset or cf.get_cfg_sec_dict(cfg, 'offset', convert='float')
    conversion = conversion or cf.get_cfg_sec_dict(cfg, 'conversion',
                                                   convert='float')
    inversion = inversion or cf.get_cfg_sec_dict(cfg, 'inversion',
                                                 convert='boolean')

    record = rc_class(datasource, traces, calibration, offset=offset,
                      conversion=conversion, inversion=inversion, name=name,
                      group=group, **record_pars)

    return record
