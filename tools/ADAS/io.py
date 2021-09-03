
import h5py
import numpy as np


def save_dict(elements, outputfile):
    """
    Save the ADAS data to an HDF5 file.
    """
    with h5py.File(outputfile, 'w') as f:
        _save_internal(elements, f)


def _save_internal(dct, f, path=''):
    """
    Internal function for saving data to HDF5.
    """
    for k in dct.keys():
        if type(dct[k]) == dict:
            o = f.create_group(k)
            _save_internal(dct[k], o, path=path+'/'+k)
        elif type(dct[k]) == float:
            f.create_dataset(k, (1,), data=dct[k])
        elif type(dct[k]) == int:
            f.create_dataset(k, (1,), data=dct[k], dtype='i8')
        elif type(dct[k]) == float:
            v = 1 if dct[k] else 0
            f.create_dataset(k, (1,), data=v, dtype='i4')
        elif type(dct[k]) == str:
            dset = f.create_dataset(k, (1,), dtype='S'+str(len(dct[k])))
            dset[0:l] = np.string_(dct[k])
        elif type(dct[k]) == list:
            f.create_dataset(k, (len(dct[k]),), data=dct[k])
        elif type(dct[k]) == np.ndarray:
            f.create_dataset(k, dct[k].shape, data=dct[k])
        else:
            raise Exception("Unrecognized data type of entry '{}/{}': {}.".format(path, k, type(dct[k])))


