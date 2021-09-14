# Unknown quantity used by the DREAMEqsys class

class DREAMEqsysUnknown:
    

    def __init__(self, id, offset, name, size, description):
        """
        Constructor.
        """
        self.id = id
        self.offset = offset
        self.name  = name
        self.size  = size
        self.description = description

        self.nr = None
        self.np = None
        self.nxi = None
        self.nZ0 = None
        self.nions = None


    def __repr__(self):
        return self.__str__()


    def __str__(self):
        """
        String representation.
        """
        return self.name


    def at(self, index):
        return (self.offset <= index and index < self.offset+self.size)


    def getNextOffset(self):
        """
        Returns the index offset for the next unknown, following this one.
        """
        return self.offset + self.size


    def getr(self, index):
        """
        Returns the radial index corresponding to the given global matrix index.
        """
        if self.nr is None:
            return None

        idx = index-self.offset
        sz  = self.size

        if self.nions is not None:
            idx -= int(self.getion(index) * (self.size/self.nions))
            sz  /= self.nions
        elif self.nZ0 is not None:
            idx -= int(self.getZ0(index) * (self.size/self.nZ0))
            sz  /= self.nZ0
            
        return int(idx / (sz / self.nr))

    
    def getp(self, index):
        """
        Returns the momentum index corresponding to the given global matrix index.
        """
        if self.np is None:
            return None

        idx = index-self.offset

        idx -= int(self.getr(index) * (self.size/self.nr))
        idx -= int(self.getxi(index) * (self.size/self.nr/self.nxi))

        return int(idx)


    def getxi(self, index):
        """
        Returns the momentum index corresponding to the given global matrix index.
        """
        if self.nxi is None:
            return None

        idx = index-self.offset

        idx -= int(self.getr(index) * (self.size/self.nr))

        return int(idx / (self.size / self.nr/self.nxi))


    def getion(self, index):
        """
        Returns the ion index corresponding to the given global matrix index.
        """
        if self.nions is None:
            return None

        idx = int((index-self.offset) / (self.size/self.nions))

        return idx


    def getZ0(self, index):
        """
        Returns the charge state index corresponding to the given global matrix index.
        """
        if self.nZ0 is None:
            return None

        idx = int((index-self.offset) / (self.size/self.nZ0))

        return idx


    def setSizes(self, nions=None, nZ0=None, nr=None, nxi=None, np=None):
        """
        Sets the size of this unknown in various dimensions.
        """
        self.nions = nions
        self.nZ0 = nZ0
        self.nr = nr
        self.nxi = nxi
        self.np = np


    def toxstring(self, index):
        """
        Convert to extended string. Includes name of unknown as
        well as (Z0,r,xi,p) indices.

        :param int index: Long format index into matrix/vector.
        """
        s = '{} ('.format(self.name)

        ii = self.getion(index)
        iz = self.getZ0(index)
        ir = self.getr(index)
        ix = self.getxi(index)
        ip = self.getp(index)

        if ii is not None: s += 'ion={}, '.format(ii)
        if iz is not None: s += 'Z0={}, '.format(iz)
        if ir is not None: s += 'r={}, '.format(ir)
        if ix is not None: s += 'xi={}, '.format(ix)
        if ip is not None: s += 'p={}, '.format(ip)

        return s[:-2] + ')'


