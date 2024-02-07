# An object representing numeric data. It is a wrapper around numpy.array and
# HDF5 data intended to allow loading data either directly into memory, or read
# it on-demand from an HDF5 file.

import h5py
import numpy as np


DATA_TYPE_ARRAY    = 1
DATA_TYPE_HDF5     = 2
DATA_TYPE_STRING   = 3
DATA_TYPE_H5STRING = 4

    
class DataObject:
    
    def __init__(self, obj):
        """
        Constructor.
        """
        self.data = obj

        if type(obj) == np.ndarray:
            self.type  = DATA_TYPE_ARRAY
            self.shape = obj.shape
            self.ndim  = obj.ndim
        elif type(obj) == h5py.Dataset:
            if obj.dtype == 'S1' or str(obj.dtype).startswith('|S') or obj.dtype == 'object':
                self.type = DATA_TYPE_H5STRING
                self.shape = (len(self[:]),)
                self.ndim  = 1
            else:
                self.type  = DATA_TYPE_HDF5
                self.shape = obj.shape
                self.ndim  = obj.ndim
        elif type(obj) == str:
            self.type  = DATA_TYPE_STRING
            self.shape = (len(obj),)
            self.ndim  = 1
        else:
            raise Exception("Unrecognized data type of object: '{}'.".format(type(obj)))

        self.size = np.prod(self.shape)


    def __eq__(self, other):
        """
        Test for equality between this object and the given
        object.
        """
        if self.type in [DATA_TYPE_ARRAY, DATA_TYPE_STRING]:
            return (self.data == other)
        else:
            return (self.data[:] == other)


    def __getitem__(self, index):
        """
        Access data.
        """
        if self.type == DATA_TYPE_H5STRING:
            # Convert string properly
            return self._getstring()[index]
        elif self.type == DATA_TYPE_HDF5:
            if type(index) == tuple and None in index:
                return self.data[:][index]
            else:
                return self.data[index]
        else:
            return self.data[index]


    def __int__(self):
        if self.type == DATA_TYPE_ARRAY:
            return int(self.data)
        else:
            return int(self[:])
    

    def __len__(self):
        """
        Return the number of elements in this dataset.
        """
        return self.size


    def __repr__(self):
        """
        String representation of this object.
        """
        if self.type == DATA_TYPE_STRING:
            return self.data.__repr__()
        elif self.type == DATA_TYPE_H5STRING:
            return self._getstring().__repr__()
        elif self.type == DATA_TYPE_ARRAY:
            return self.data.__repr__()
        else:
            return self.data[:].__repr__()


    def __str__(self):
        """
        Convert this data object into a string.
        """
        if self.type == DATA_TYPE_ARRAY:
            return str(self.data)
        else:
            return str(self[:])


    def _getstring(self):
        """
        Convert an H5STRING to a Python string.
        """
        if self.type == DATA_TYPE_H5STRING:
            # Convert string properly
            if (self.data.dtype == 'S1') or (str(self.data.dtype).startswith('|S')):  # Regular strings
                return self.data[:].tostring().decode('utf-8')
            elif self.data.dtype == 'object':  # New strings
                if self.data.shape == ():
                    return self.data[()].decode()
                elif type(self.data[:][0]) == str:
                    return self.data[:][0]
                else:
                    return self.data[:][0].decode()
        elif self.type == DATA_TYPE_STRING:
            return self.data
        else:
            raise Exception("The '_getstring()' method is only intended for HDF5 strings.")


    def __array__(self, dtype=None):
        """
        Convert to numpy array.
        """
        return np.asarray(self[:], dtype=dtype)


    """
    ARITHMETIC OPERATIONS
    """
    def _lop(self, op, f):
        """
        Binary operator where self is the left operand.
        """
        if np.isscalar(op):
            return f(self.data[:], op)
        else:
            return f(self.data[:], op[:])


    def __add__(self, other): return self._lop(other, lambda o1, o2 : o1+o2)


    def __sub__(self, other): return self._lop(other, lambda o1, o2 : o1-o2)


    def __mul__(self, other): return self._lop(other, lambda o1, o2 : o1*o2)


    def __truediv__(self, other): return self._lop(other, lambda o1, o2 : o1.__truediv__(o2))


    def __pow__(self, other): return self._lop(other, lambda o1, o2 : o1**o2)


    def __and__(self, other): return self._lop(other, lambda o1, o2 : o1 and o2)


    def __or__(self, other): return self._lop(other, lambda o1, o2 : o1 or p2)


    def __neg__(self): return -self.data[:]


    def __pos__(self): return +self.data[:]


    def __abs__(self): return abs(self.data[:])


    def __invert__(self): return ~self.data[:]


