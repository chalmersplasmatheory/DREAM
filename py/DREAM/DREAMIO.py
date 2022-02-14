#
# Common file interface.
# #######################

import h5py
import numpy as np
from packaging import version
from pathlib import Path
import os

from . DataObject import DataObject

# Try to import paramiko for SSH support (optional)
SSHSUPPORT = False
try:
    import paramiko
    import re           # Regular expression matcher
    SSHSUPPORT = True
except ModuleNotFoundError: pass


def LoadHDF5AsDict(filename, path='', returnhandle=False, returnsize=False, lazy=True):
    """
    Loads the given HDF5 file as a dict.

    :param str filename:      Name of HDF5 file to load.
    :param str path:          Path to subset of HDF5 file to load.
    :param bool returnhandle: If ``True``, also returns the HDF5 file handled
                              used to open the file (always ``None`` if ``lazy==False``).
    :param bool returnsize:   If ``True``, also returns the file size.
    :param bool lazy:         If ``True``, loads data as ``DataObject``, which
                              allows lazy (on-demand) reading of data. The
                              default is to immediately load all data from the
                              file.
    """
    global SSHSUPPORT
    data = None
    size = 0
    f = None

    user, host, port, rpath = None, None, 22, None
    if SSHSUPPORT:
        m1 = re.search('(\w+://)(.+@)*([\w\-\_\d\.]+)(:[\d]+){0,1}/*(.*)', filename)
        m2 = re.search('(.+@)*([\w\-\_\d\.]+):(.*)', filename)

        if m1 is not None:
            user = m1.group(2)
            host = m1.group(3)
            prtt = m1.group(4)
            rpath = m1.group(5)

            # Remove '@' in username (if given)
            if user is not None:
                user = user[:-1]

            if prtt is not None:
                port = int(prtt[1:])
            else:
                rpath = rpath[1:]

        elif m2 is not None:
            user = m2.group(1)
            host = m2.group(2)
            rpath = m2.group(3)

            # Remove '@' in username (if given)
            if user is not None:
                user = user[:-1]

    if host is not None and rpath is not None:
        client = paramiko.SSHClient()
        client.load_system_host_keys()

        config = paramiko.config.SSHConfig()
        try:
            config.parse(open("{}/.ssh/config".format(str(Path.home()))))

            conf = config.lookup(host)
            host = conf['hostname']

            if 'user' in conf:
                user = conf['user']
            if 'port' in conf:
                port = conf['port']

        except: pass

        client.connect(host, port=port, username=user)

        # Open SFTP stream
        sftp = client.open_sftp()
        size = sftp.stat(rpath).st_size
        
        if lazy:
            fo = sftp.open(rpath, 'r')
            f  = h5py.File(fo, 'r')
            data = h52dict(f, path, lazy=True)
            # Close neither connection nor HDF5 file to allow lazy reading...
        else:
            with sftp.open(rpath, 'r') as fo:
                with h5py.File(fo, 'r') as f:
                    data = h52dict(f, path, lazy=False)
            client.close()
    else:
        size = os.path.getsize(filename)
        if lazy:
            f = h5py.File(filename, 'r')
            data = h52dict(f, path, lazy=True)
            # Don't close HDF5 file to allow lazy reading...
        else:
            with h5py.File(filename, 'r') as f:
                data = h52dict(f, path, lazy=False)

    ret = (data, )
    if returnhandle:
        if lazy:
            ret += (f,)
        else:
            ret += (None,)
    if returnsize: ret += (size,)

    if len(ret) > 1:
        return ret
    else:
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
        d = data[key]
        if type(d) == DataObject:
            d = d[:]

        if type(d) == dict:
            o = f.create_group(key)
            dict2h5(o, d, path=path+'/'+key)
        elif type(d) == float:
            f.create_dataset(key, (1,), data=d)
        elif type(d) == int:
            f.create_dataset(key, (1,), data=d, dtype='i8')
        elif type(d) == bool:
            v = 1 if d else 0
            f.create_dataset(key, (1,), data=v, dtype='i4')
        elif type(d) == str:
            l = len(d)

            # From h5py version 2.10.0 an on there is support for storing
            # UTF-8 strings. To allow this, but still remain backwards
            # compatible, we choose how to store strings depending on the
            # version of h5py.
            if version.parse(h5py.version.version) >= version.parse('2.10.0'):
                dt = h5py.string_dtype()
                dset = f.create_dataset(key, (1,), dtype=dt)
                dset[0:l] = d
            else:   # h5py < 2.10.0
                dset = f.create_dataset(key, (1,), dtype='S'+str(l))
                dset[0:l] = np.string_(d)
        elif type(d) == list:
            f.create_dataset(key, (len(d),), data=d)
        elif type(d) == np.ndarray:
            f.create_dataset(key, d.shape, data=d)
        else:
            raise DREAMIOException("Unrecognized data type of entry '{}/{}': {}.".format(path, key, type(d)))


def h52dict(f, path='', lazy=False):
    """
    Loads data from the given HDF5 file handle 'f'.

    f:    HDF5 file handle to use for reading.
    path: Path in HDF5 file to read data from.
    lazy: Load data lazily, i.e. return a DataObject rather than
          the actual data.
    """
    d = {}
    for key in f.keys():
        if type(f[key]) == h5py.Group:
            d[key] = h52dict(f[key], path=path+'/'+key, lazy=lazy)
        elif type(f[key]) == h5py.Dataset:
            if lazy:
                d[key] = DataObject(f[key])
            else:
                d[key] = getData(f, key)

            # Get attributes (cannot be loaded lazily)
            if len(f[key].attrs) > 0:
                n = key+'@@'
                if n not in d:
                    d[n] = {}

                for a in f[key].attrs:
                    d[key+'@@'][a] = getData(f[key].attrs, a)
        else:
            raise DREAMIOException("Unrecognized HDF5 data structure for key: '{}'.".format(path+'/'+key))

    return d


def getData(f, key):
    """
    Returns data from an h5py.File object, correctly transforming
    it (in case it is a string for example).
    """
    if (f[key].dtype == 'S1') or (str(f[key].dtype).startswith('|S')):  # Regular strings
        return f[key][:].tostring().decode('utf-8')
    elif f[key].dtype == 'object':  # New strings
        if f[key].shape == ():
            return f[key][()].decode()
        else:
            return f[key][:][0].decode()
    else:
        return f[key][:]


class DREAMIOException(Exception):
    def __init__(self, msg):
        super(Exception, self).__init__(msg)


