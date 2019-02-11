# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:04:15 2016

@author: Tobias Jachowski
"""
import numpy as np
#import os
import struct

'''
def get_labview_bin_data_shape(filename):
    print('Getting data shape from:')
    print('  \'%s\'' % filename)
    filesize = os.path.getsize(filename) / 8  # size in bytes
    chunk_shape = np.fromfile(filename, count=2, dtype='>i4')
    # size of buffer, which is read in every loop iteration
    buffer_size = chunk_shape[0]
    # number of channels read in
    number_channels = chunk_shape[1]
    # chunk size
    chunk_size = 1 + buffer_size * number_channels
    # chunk_shape of single chunks is not changed during a measurement with the
    # ForceClamp -> same chunk_size for all chunks
    # number of chunks in file
    number_chunks = int(filesize/(chunk_size))
    number_samples = buffer_size * number_chunks

    shape_file = (number_chunks, chunk_size)
    shape_record = (number_samples, number_channels)

    return shape_file, shape_record


def read_labview_bin_data(filename, dtype='>d'):
    # big-endian double-precision floating-point number
    shape_file, shape_record = get_labview_bin_data_shape(filename)

    # 1. read in data
    # 2. reshape to two dimensions: chunks and bytes(/samples)
    # 3. remove byte corresponding to the chunk shape information
    # 4. reshape data to: number of samples and number of channels
    count = shape_file[0] * shape_file[1]
    # workaround, if file filesize % chunk_size is not an integer
    # it leaves out additional numbers

    print("Reading data from:\n  '%s'" % filename)
    data = np.fromfile(
        filename,
        dtype=dtype,
        count=count).reshape(shape_file)[:, 1:].reshape(shape_record)
    # data = np.loadtxt(filename, comments='%', delimiter=' ')

    # Workaround, if last column is 0.0
    if np.all(data[:, -1] == 0.0):
        # kick out last column
        data = data[:, :-1]

    # data = np.memmap(filename,
    #                  dtype=dtype,
    #                  mode='r', shape=(l, k))[:,1:].reshape([l*i, j])
    #
    # fname_o = ''
    # convert_bin_data(filename, fname_o, k, 1)
    # data = np.memmap(fname_o, dtype=dtype, mode='r', shape=(l*i, j))
    # data = np.fromfile(fname_o, dtype=dtype).reshape([l*i, j])
    return data
'''


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


def __OLD_read_labview_bin_data(filename, dtype='>d'):
    rows = 0
    columns = None
    chunks = []

    print('Getting chunk info from:')
    print('  \'%s\'' % filename)
    for chunk in chunk_info(filename):
        c = chunk[2]
        if columns == c or columns is None:
            columns = c
        else:
            print("Number of columns in chunks of file differ from each other!")
        rows += chunk[1]
        chunks.append(chunk)

    print('Reading data chunks from:')
    print('  \'%s\'' % filename)
    data = np.empty(rows * columns)
    i = 0
    for _data in chunk_data(filename, chunks, dtype=dtype):
        length = len(_data)
        data[i:i + length] = _data
        i += length

    data = data.reshape(rows, columns)

    # Workaround, if last column is 0.0
    if np.all(data[:, -1] == 0.0):
        # kick out last column
        data = data[:, :-1]

    return data


def write_labview_bin_data(filename, array, dtype='d', endianness='>'):

    i = struct.pack('>i', array.shape[0])
    j = struct.pack('>i', array.shape[1])
    fmt = endianness + dtype * array.size
    data = struct.pack(fmt, *array.flatten())

    with open(filename, 'wb') as f:
        f.write(i)
        f.write(j)
        f.write(data)

'''
def convert_bin_data(fname_i, fname_o, chunk_size, chunk_shape_size=1):
    """
    Remove chunk size information of labview binary files.
    """
    fi = open(fname_i, 'rb')
    fo = open(fname_o, 'wb')

    while True:
        buf = fi.read(chunk_size)
        if len(buf) == 0:
            break
        fo.write(buf[chunk_shape_size:])

        fi.close()
        fo.close()
'''