# DREAM Equation System parser
#
# This classes parses the equation system header written to stdout by DREAM when
# a simulation is started. The class then allows the user to convert from matrix
# and vector indices to DREAM unknown quantity names and (Z0,r,xi,p) indices.
#
#################

import numpy as np
from .DREAMEqsysUnknown import DREAMEqsysUnknown


class DREAMEqsys:
    """
    This class loads the equation system stdout output from a DREAM
    simulation, parses it and allows you to map unknowns by name to
    matrix indices.
    """

    def __init__(self, filename, nions=None):
        """
        Constructor.

        :param nions: Number of ion species in simulation (not required if 'N_i' and/or 'W_i' is present in the equation system).
        """
        self.unknowns = []
        self.nr = None
        self.nions = nions
        self.nZ0 = None
        self.hot_npnxi = None
        self.hot_np = None
        self.hot_nxi = None
        self.re_npnxi = None
        self.re_np = None
        self.re_nxi = None

        with open(filename, 'r') as f:
            for line in f:
                l = [s for s in line.strip().split() if s]

                id = int(l[0])
                name  = l[1]
                size  = int(l[2])
                descr = ' '.join(l[3:])

                self.add(id=id, name=name, size=size, description=descr)

        self.initialize_sizes()


    def __contains__(self, name):
        """
        Check if this equation system contains the named unknown.
        """
        for u in self.unknowns:
            if u.name == name:
                return True

        return False


    def __getitem__(self, index):
        """
        Returns the unknown with the given name or at the specified
        matrix index.
        """
        if type(index) == str:
            return self.get(index)
        else:
            return self.at(index)


    def at(self, index):
        """
        Returns the name of the unknown at the specified index.
        """
        for u in self.unknowns:
            if u.at(index):
                return u

        return None


    def add(self, id, name, size, description):
        """
        Add an unknown to the equation system.
        """
        if len(self.unknowns) > 0:
            offset = self.unknowns[-1].getNextOffset()
        else:
            offset = 0

        self.unknowns.append(DREAMEqsysUnknown(id=id, offset=offset, name=name, size=size, description=description))


    def get(self, name):
        """
        Return the unknown with the given name.
        """
        for u in self.unknowns:
            if u.name == name:
                return u

        return None


    def getnames(self):
        """
        Return the name of each unknown in a list.
        """
        return [u.name for u in self.unknowns]


    def getr(self, index):
        """
        Returns the radial index corresponding to the given global matrix index.
        """
        return self.at(index).getr(index)

    
    def getp(self, index):
        """
        Returns the momentum index corresponding to the given global matrix index.
        """
        return self.at(index).getp(index)


    def getxi(self, index):
        """
        Returns the pitch index corresponding to the given global matrix index.
        """
        return self.at(index).getxi(index)


    def getion(self, index):
        """
        Returns the ion index corresponding to the given global matrix index.
        """
        return self.at(index).getion(index)


    def getZ0(self, index):
        """
        Returns the charge state index corresponding to the given global matrix index.
        """
        return self.at(index).getZ0(index)


    def getoffsets(self):
        """
        Return the offset of each unknown in a list.
        """
        return [u.offset for u in self.unknowns]


    def initialize_sizes(self):
        """
        Determine grid resolution.
        """
        # Determine NR
        fluid = ['n_cold', 'n_hot', 'n_re', 'E_field', 'T_cold']
        for f in fluid:
            if f in self:
                self.nr = self[f].size
                break

        # Number of charge states
        if 'n_i' in self:
            self.nZ0 = self['n_i'].size / self.nr

        # Number of ion species
        if 'N_i' in self:
            self.nions = self['N_i'].size / self.nr
        elif 'W_i' in self:
            self.nions = self['W_i'].size / self.nr

        # Hot np*nxi
        if 'f_hot' in self:
            self.hot_npnxi = self['f_hot'].size / self.nr

        # Runaway np*nxi
        if 'f_re' in self:
            self.re_npnxi = self['f_re'].size / self.nr

        self.update_unknown_sizes()


    def setHot(self, np=None, nxi=None):
        """
        Set hot grid resolution. Usually, only one of the parameters needs to be
        specified; the other one can be deduced from the size of 'f_hot'.
        """
        if np is None and nxi is None:
            return

        if np is not None:
            self.hot_np = np
            self.hot_nxi = self.hot_npnxi / np

        if nxi is not None:
            self.hot_nxi = nxi
            
            if self.hot_np is None:
                self.hot_np = self.hot_npnxi / nxi

        if self.hot_np*self.hot_nxi != self.hot_npnxi:
            raise Exception("Invalid size of f_hot specified. np*nxi should be {}.".format(self.hot_npnxi))

        self.update_unknown_sizes()


    def setRE(self, np=None, nxi=None):
        """
        Set hot grid resolution. Usually, only one of the parameters needs to be
        specified; the other one can be deduced from the size of 'f_re'.
        """
        if np is None and nxi is None:
            return

        if np is not None:
            self.re_np = np
            self.re_nxi = self.re_npnxi / np

        if nxi is not None:
            self.re_nxi = nxi
            
            if self.re_np is None:
                self.re_np = self.re_npnxi / nxi

        if self.re_np*self.re_nxi != self.re_npnxi:
            raise Exception("Invalid size of f_re specified. np*nxi should be {}.".format(self.re_npnxi))

        self.update_unknown_sizes()


    def update_unknown_sizes(self):
        """
        Updates the size variables of all unknowns.
        """
        for u in self.unknowns:
            if u.name == 'f_hot':
                u.setSizes(nr=self.nr, nxi=self.hot_nxi, np=self.hot_np)
            elif u.name == 'f_re':
                u.setSizes(nr=self.nr, nxi=self.re_nxi, np=self.re_np)
            elif u.name == 'n_i':
                u.setSizes(nZ0=self.nZ0, nr=self.nr)
            elif u.size == self.nr:     # Fluid quantity
                u.setSizes(nr=self.nr)
            elif u.size == self.nions*self.nr:      # N_i or W_i
                u.setSizes(nions=self.nions, nr=self.nr)


    def __str__(self):
        """
        String representation.
        """
        s = "ID   NAME                START   DESCRIPTION\n"
        for u in self.unknowns:
            s += "{:3d}  {:15s} {:9d}   {}\n".format(u.id, u.name, u.offset, u.description)

        return s



