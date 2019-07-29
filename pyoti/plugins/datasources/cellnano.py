# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:37:58 2016

@author: Tobias Jachowski
"""
import numpy as np
import struct
import os

from pyoti.data.datasource import DataSource


class CNLabViewBinData(DataSource):
    def __init__(self, filename, directory=None, parafile=None, datext='.bin',
                 parext='_para.dat', **kwargs):
        """
        parext : str, optional
            The extension of the parameter file ('_para.dat', default)
        """
        super().__init__(filename=filename, directory=directory, **kwargs)

        if parafile is None:
            parafile = self.absfile.replace(datext, parext)
        self._parafile_orig = parafile

        self.samplingrate = get_samplingrate_from_parameter_file(parafile)

        self.name = ("LabVIEW bin data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (self.absfile_orig,
                                                  self.samplingrate)

    def as_array(self):
        filename = self.absfile
        data = read_labview_bin_data(filename)
        return data

    @property
    def parafile_orig(self):
        return self._parafile_orig


class CNLabViewTxtData(DataSource):
    def __init__(self, filename, directory=None, samplingrate=1000.0,
                 **kwargs):
        """
        parext : str, optional
            The extension of the parameter file ('_para.dat', default)
        """
        super().__init__(filename=filename, directory=directory, **kwargs)
        self.samplingrate = samplingrate

        self.name = ("LabVIEW txt data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (self.absfile_orig,
                                                  self.samplingrate)

    def as_array(self):
        filename = self.absfile
        data = np.loadtxt(filename, skiprows=5)
        return data


class CNLabViewQPDviData(CNLabViewTxtData):
    def __init__(self, filename, directory=None, **kwargs):
        super().__init__(filename, directory=directory, **kwargs)
        self.samplingrate = np.genfromtxt(self.absfile, skip_header=2,
                                          max_rows=1)[0]

        self.name = ("LabVIEW QPD.vi data originally loaded from \n"
                     "    %s with \n"
                     "    samplingrate %s Hz") % (os.path.join(self.directory,
                                                               self.filename),
                                                  self.samplingrate)


def get_samplingrate_from_parameter_file(parafile):
    para = np.loadtxt(parafile, comments='%', delimiter='\t')
    scanrate = para[0]  # float
    decimating = para[2]  # float
    samplingrate = scanrate / decimating  # float
    return samplingrate


def chunk_info(filename):
    """
    Read all chunk arrays with the information (number of rows and columns) of
    the data, following each chunk.
    """
    with open(filename, "rb") as f:
        pos = 0
        while True:
            byte = f.read(8)
            pos += 1
            if byte:
                chunk_shape = struct.unpack('>2i', byte)
                if chunk_shape[0] * chunk_shape[1] > 0:
                    # go to byte offset relative to current strem position
                    f.seek(8 * chunk_shape[0] * chunk_shape[1], 1)
                    yield(pos, *chunk_shape)
                    pos += chunk_shape[0] * chunk_shape[1]
            else:
                break


def chunk_data(filename, chunks, dtype='>d'):
    with open(filename, "rb") as f:
        for chunk in chunks:
            # go to byte offset where chunk starts
            f.seek(chunk[0] * 8)
            # read the data
            data_bin = f.read(chunk[1] * chunk[2] * 8)
            data = np.fromstring(data_bin, dtype=dtype)
            yield data


def read_labview_bin_data(filename, dtype='>d', start_row_idx=0,
                          number_of_rows=-1):
    """
    Parameters
    ----------
    filename : str
        Path of the binary file to read
    dtype : str
        Type of double
    start_row_idx : int
        Index of the first datapoint (of all traces) to read
    number_or_rows : int
        Number of datapoints (of all traces) to read. Defaults to number of
        datapoints of the binary file - `start_row_idx`
    """
    stop_row_idx = start_row_idx + number_of_rows  # index of row to stop read
    rows_to_read = 0  # rows to read after getting information of chunks
    columns = None  # number of columns (i.e. traces) in binary file
    chunks = []  # chunks with the index information of data to be read
    chunk_row_start = 0  # running index of first row of chunk
    chunk_row_stop = 0  # running stop index of last row of chunk

    if number_of_rows == 0:
        return np.empty((0, 0))

    print('Getting chunk info from:')
    print('  \'%s\'' % filename)
    for chunk in chunk_info(filename):
        # Set the running indices to the new position (row)
        chunk_row_start = chunk_row_stop
        chunk_row_stop += chunk[1]
        # Check if information about number of columns changed from one to the
        # next chunk
        if columns == chunk[2] or columns is None:
            columns = chunk[2]
        else:
            print("Number of columns in chunks of file differ from each other!")
        # Check if chunk is (partly) contained within the requested data index
        if (chunk_row_stop > start_row_idx
                and (number_of_rows < 0 or chunk_row_start < stop_row_idx)):
            # Check read position and number of rows of first chunk
            if len(chunks) is 0:
                shift = start_row_idx - chunk_row_start
                chunk = (chunk[0] + shift * chunk[2],
                         chunk[1] - shift,
                         chunk[2])
            chunks.append(chunk)
            rows_to_read += chunk[1]
        # Check stop position and number of rows of last chunk
        if (number_of_rows > 0 and chunk_row_stop >= stop_row_idx):
            shift = chunk_row_stop - stop_row_idx
            chunk = (chunk[0], chunk[1] - shift, chunk[2])
            chunks[-1] = chunk
            rows_to_read -= shift
            break

    if rows_to_read == 0:
        return np.empty((0, 0))

    print('Reading data chunks from:')
    print('  \'%s\'' % filename)
    data = np.empty(rows_to_read * columns)
    i = 0
    for _data in chunk_data(filename, chunks, dtype=dtype):
        length = len(_data)
        data[i:i + length] = _data
        i += length

    data = data.reshape(rows_to_read, columns)

    # Workaround, if last column is 0.0
    if np.all(data[:, -1] == 0.0):
        # kick out last column
        data = data[:, :-1]

    return data
