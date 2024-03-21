# Class for diffing two settings files

import h5py
import numpy as np


class SettingsDifference:
    

    def __init__(self, path, reason, val1, val2, type1=None, type2=None, shape1=None, shape2=None):
        """
        Constructor.
        """
        self.path = path
        self.reason = reason
        self.val1 = val1
        self.val2 = val2
        self.type1 = type1
        self.type2 = type2
        self.shape1 = shape1
        self.shape2 = shape2


class SettingsDiff:
    

    def __init__(self):
        """
        Constructor.
        """
        pass


    def diff(self, file1, file2):
        """
        Make a diff between the two files.
        """
        self.file1 = file1
        self.path1 = '/'
        self.file2 = file2
        self.path2 = '/'

        if ':' in file1:
            ff = file1.split(':')
            self.file1 = ff[0]
            self.path1 = ff[-1]
        if ':' in file2:
            ff = file2.split(':')
            self.file2 = ff[0]
            self.path2 = ff[-1]

        self.similarities = 0
        el = []
        with h5py.File(self.file1) as f1, h5py.File(self.file2) as f2:
            el += self._diff(f1[self.path1], f2[self.path2], path='')

        return el


    def _diff(self, f1, f2, path, ignore_nonexisting=False):
        """
        Diff the two HDF5 datasets.
        """
        if type(f1) != type(f2):
            return [SettingsDifference(path, 'Different types', f'Has type {type(f1)}', f'Has type {type(f2)}')]
        elif type(f1) == h5py.Group:
            diffs = []
            k1 = list(f1.keys())
            k2 = list(f2.keys())
            if len(k1) != len(k2):
                for k in k2:
                    if k not in k1:
                        if not ignore_nonexisting:
                            diffs.append(SettingsDifference(f'{path}/{k}', 'Non-existant groups', 'NOT EXISTS', 'EXISTS'))
            for k in k1:
                if k not in k2:
                    if not ignore_nonexisting:
                        diffs.append(SettingsDifference(f'{path}/{k}', 'Non-existant groups', 'EXISTS', 'NOT EXISTS'))
                else:
                    diffs += self._diff(f1[k], f2[k], path=f'{path}/{k}', ignore_nonexisting=ignore_nonexisting)

            return diffs
        else:
            if f1.shape != f2.shape:
                return [
                    SettingsDifference(
                        path, 'Different shapes',
                        f'{f1[:]}', f'{f2[:]}',
                        type1=f'{type(f1)}', type2=f'{type(f2)}',
                        shape1=f1.shape, shape2=f2.shape
                    )
                ]
            elif np.all(f1[:] == f2[:]):
                self.similarities += 1
                return []
            else:
                return [
                    SettingsDifference(
                        path, 'Different values',
                        f'{f1[:]}', f'{f2[:]}',
                        type1=f'{type(f1)}', type2=f'{type(f2)}',
                        shape1=f1.shape, shape2=f2.shape
                    )
                ]


