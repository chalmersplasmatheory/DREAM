#
# Common file interface.
# #######################

import h5py
import numpy as np


def LoadHDF5AsDict(filename):
    """
    Loads the given HDF5 file as a dict.

    filename: Name of HDF5 file to load.
    """
    data = None
    with h5py.File(filename, 'r') as f:
        data = h52dict(f)

    return data


def SaveDictAsHDF5(filename, data):
    """
    Stores the given dict into an HDF5 file.

    filename: Name of HDF5 file to create and write to.
    data:     Dictionary containing the data to write.
    """
    with h5py.File(filename, 'w') as f:
        dict2h5(f, data)


def dict2h5(f, data, path=''):
    """
    Writes the given data on the HDF5 file handle 'f'.

    f:    HDF5 file handle to write on.
    data: Data to write.
    path: Path in HDF5 file to write data to.
    """

    for key in data:
        if type(data[key]) == dict:
            o = f.create_group(key)
            dict2h5(o, data[key], path=path+'/'+key)
        elif type(data[key]) == float:
            f.create_dataset(key, (1,), data=data[key])
        elif type(data[key]) == int:
            f.create_dataset(key, (1,), data=data[key], dtype='i8')
        elif type(data[key]) == bool:
            v = 1 if data[key] else 0
            f.create_dataset(key, (1,), data=v, dtype='i4')
        elif type(data[key]) == str:
            l = len(data[key])
            dset = f.create_dataset(key, (1,), dtype='S'+str(l))
            dset[0:l] = np.string_(data[key])
        elif type(data[key]) == list:
            f.create_dataset(key, (len(data[key]),), data=data[key])
        elif type(data[key]) == np.ndarray:
            f.create_dataset(key, data[key].shape, data=data[key])
        else:
            raise DREAMIOException("Unrecognized data type of entry '{}/{}'.".format(path, key))


def h52dict(f, path=''):
    """
    Loads data from the given HDF5 file handle 'f'.

    f:    HDF5 file handle to use for reading.
    path: Path in HDF5 file to read data from.
    """
    d = {}
    for key in f.keys():
        if type(f[key]) == h5py.Group:
            d[key] = h52dict(f[key], path=path+'/'+key)
        elif type(f[key]) == h5py.Dataset:
            d[key] = f[key][:]
        else:
            raise DREAMIOException("Unrecognized HDF5 data structure for key: '{}'.".format(path+'/'+key))

    return d


class DREAMIOException(Exception):
    def __init__(self, msg):
        super(Exception, self).__init__(msg)


