# Settings object for the hot-tail grid.

import numpy as np
from DREAM.DREAMException import DREAMException


from DREAM.Settings.PGrid import PGrid
from DREAM.Settings.XiGrid import XiGrid


TYPE_PXI = 1
TYPE_PPARPPERP = 2


class MomentumGrid:

    
    def __init__(self, name, enabled=True, ttype=TYPE_PXI, np=0, nxi=0, pmax=None):
        """
        Constructor.

        name:    Name, indicating what type of grid this is (hot-tail or runaway).
        enabled: If 'True', enables the hot-tail grid in the simulation.
        ttype:   Type of momentum grid (p/xi or ppar/pperp).
        np:      Number of momentum grid points.
        nxi:     Number of pitch grid points.
        pmax:    Maximum momentum on grid.
        
        """
        self.name = name

        self.set(enabled=enabled, ttype=ttype, np=np, nxi=nxi, pmax=pmax)


    def __contains__(self, item):
        return (item in self.todict(False))


    def __getitem__(self, index):
        return self.todict(False)[index]


    def set(self, enabled=True, ttype=1, np=100, nxi=1, pmax=None):
        """
        Set all settings for this hot-tail grid.
        """
        self.enabled = enabled

        self.type = ttype
        if self.type == TYPE_PXI:
            self.pgrid = PGrid(self.name, np=np, pmax=pmax)
            self.xigrid = XiGrid(self.name, nxi=nxi)
        elif self.type == TYPE_PPARPPERP:
            raise DREAMException("No support implemented yet for 'ppar/pperp' grids.")
        else:
            raise DREAMException("Unrecognized momentum grid type specified: {}.".format(ttype))

    ##################
    # SETTERS
    ##################
    def setEnabled(self, enabled=True):
        self.enabled = (enabled == True)


    def setNp(self, np):
        if np <= 0:
            raise DREAMException("{}: Invalid value assigned to 'np': {}. Must be > 0.".format(self.name, np))
        elif np == 1:
            print("WARNING: {}: np = 1. Consider disabling the hot-tail grid altogether.".format(self.name))

        self.pgrid.setNp(np)


    def setNxi(self, nxi):
        if nxi <= 0:
            raise DREAMException("{}: Invalid value assigned to 'nxi': {}. Must be > 0.".format(self.name, nxi))

        self.xigrid.setNxi(nxi)


    def setPmax(self, pmax):
        if pmax <= 0:
            raise DREAMException("{}: Invalid value assigned to 'pmax': {}. Must be > 0.".format(self.name, pmax))

        self.pgrid.setPmax(pmax)

    def setBiuniformGrid(self, psep=None, npsep=None, npsep_frac=None, xisep=None, nxisep=None, nxisep_frac=None, thetasep=None, nthetasep=None, nthetasep_frac=None):
        """
        Set a two-region momentum grid. The lower part (0 < p < psep) has 'npsep' 
        number of grid points, while the upper region has 'np-npsep' number of
        grid points. This makes it possible to, for example, have denser
        grid near the bulk and a coarser grid in the runaway tail. 
        Similarly, the lower part (-1 < xi < xisep) of the pitch grid has 'nxisep'
        number of grid points, etc.

        :param float psep:        	 Momentum value separating the two sections.
        :param int npsep:         	 Number of grid points on the lower grid section.
        :param float npsep_frac:  	 If ``npsep`` is ``None``, gives the fraction of grid 
                                  	 points to put in the lower region. Otherwise, not used.
        :param float xisep:       	 Pitch value separating the two sections.
        :param int nxisep:        	 Number of grid points on the lower grid section.
        :param float nxisep_frac: 	 If ``nxisep`` is ``None``, gives the fraction of grid
                                  	 points to put in the lower region. Otherwise not used.
        :param float thetasep:    	 Theta value separating the two sections.
        :param int nthetasep:        Number of grid points on the lower grid section.
        :param float nthetasep_frac: If ``nthetasep`` is ``None``, gives the fraction of grid
                                  	 points to put in the lower region. Otherwise not used.
        """
        if psep is not None:
            self.pgrid.setBiuniform(psep=psep,npsep=npsep,npsep_frac=npsep_frac)
        if xisep is not None:
            self.xigrid.setBiuniform(xisep=xisep,nxisep=nxisep,nxisep_frac=nxisep_frac)
        elif thetasep is not None:
            self.xigrid.setBiuniform(thetasep=thetasep,nthetasep=nthetasep,nthetasep_frac=nthetasep_frac)
            
    def setCustomGrid(self, p_f=None, xi_f=None):
        """
        Sets p and xi to an arbitrary meshgrid, where the flux grid points 
        of the respective grids are contained in the lists 'p_f' and 'xi_f'.
        The resolution parameters 'np' and 'nxi' are overwritten by the sizes
        of the two input grids.

        :param float p_f:       Increasing list of momentum flux grid points, with p_f[0]=0.
        :param float xi_f:      Increasing list of pitch grid points with xi_f[0] = -1 and xi_f[-1] = 1
        """
        if p_f is not None:
            self.pgrid.setCustomGridPoints(p_f=p_f)
        if xi_f is not None:
            self.xigrid.setCustomGridPoints(xi_f=xi_f)

    def setTrappedPassingBoundaryLayerGrid(self, xi0Trapped=None, dxiMax=2, NxiPass=1, NxiTrap=1, boundaryLayerWidth=1e-3):
        """
        Designs a custom pitch grid which places tight grid cells 
        straddling each trapped-passing boundary, so that these
        singular points are well resolved.

        :param float xi0Trapped:            List of trapped-passing boundary pitches 
                                            (contained in do.grid.xi0TrappedBoundary of an output object)
        :param float dxiMax:                Maximum allowed grid spacing dxi: if this spacing is exceeded after
                                            placing points on the trapped-passing boundary, will fill in the gaps
                                            by adding one or more grid points (uniformly) in the gaps
        :param int NxiPass:                 Number of grid cells in pitch to be contained 
                                            between the largest xi0Trapped and xi0=+/-1
        :param int NxiTrap:                 Number of grid cells in pitch to be contained 
                                            between (+/-) the minimum xi0Trapped and xi0=0
        :param float boundaryLayerWidth:    Width in pitch of the grid cell containing each
                                            trapped-passing boundary, typically << 1
        """
        if xi0Trapped is None:
            self.xigrid.setTrappedPassingBoundaryLayerGrid(dxiMax=dxiMax, NxiPass=NxiPass, NxiTrap=NxiTrap, boundaryLayerWidth=boundaryLayerWidth)
        else:
            # Old method, which requires prior knowledge of 'xi0Trapped'
            if type(xi0Trapped)==list:
                xi0Trapped = np.array(xi0Trapped)
            xi0Trapped = np.sort(xi0Trapped)

            x0 = np.linspace( xi0Trapped.max(), 1, NxiPass+1 )[1:]
            xi_f = np.array([-x0,x0])
            if NxiTrap > 1 and xi0Trapped.min() > 0:
                xTrap = np.linspace( 0, xi0Trapped.min(), NxiTrap+1 )[1:-1]
                xi_f = np.append( xi_f, [-xTrap, xTrap] )
            xi_f = np.append(xi_f, 0)

            for i in range(xi0Trapped.size):
                if xi0Trapped[i]>0:
                    xiAdd1 = xi0Trapped[i] + 0.5*boundaryLayerWidth
                    xiAdd2 = xi0Trapped[i] - 0.5*boundaryLayerWidth
                    xi_f = np.append(xi_f, [-xiAdd1, xiAdd1, -xiAdd2, xiAdd2])

            # Add additional point if the resulting spacing exceeds the desired dxiMax
            xi_f = np.sort(xi_f)
            for i in range(xi_f.size-1):
                dxi = xi_f[i+1] - xi_f[i]
                N = int(np.floor( dxi/dxiMax ))
                if N>=1:
                    xi_f = np.append(xi_f, np.linspace(xi_f[i], xi_f[i+1],N+2)[1:-1])

            self.setCustomGrid(xi_f=np.sort(xi_f))

    def setXiType(self,ttype):
        """
        Set type of xi grid. 
        """
        self.xigrid.setType(ttype)
	
    def fromdict(self, name, data):
        """
        Loads a momentum grid from the specified dictionary.
        """
        self.name    = name
        self.enabled = data['enabled']
        self.type    = data['type']

        if self.enabled:
            if self.type == TYPE_PXI:
                self.pgrid = PGrid(name, data=data)
                self.xigrid = XiGrid(name, data=data)
            elif self.type == TYPE_PPARPPERP:
                raise DREAMException("No support implemented yet for loading 'ppar/pperp' grids.")
            else:
                raise DREAMException("Unrecognized momentum grid type specified: {}.".format(self.type))
        else:
            # Set default grid
            self.set(enabled=False, ttype=self.type)

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this MomentumGrid object.
        """
        if verify:
            self.verifySettings()

        data = {
            'enabled': self.enabled,
            'type':    self.type
        }

        if not self.enabled: return data

        if self.type == TYPE_PXI:
            data = {**data, **(self.pgrid.todict()), **(self.xigrid.todict())}
        elif self.type == TYPE_PPARPPERP:
            raise DREAMException("No support implemented yet for saving 'ppar/pperp' grids.")
        else:
            raise DREAMException("Unrecognized momentum grid type specified: {}.".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that all (mandatory) settings are set and consistent.
        """
        if self.type == TYPE_PXI:
            if self.enabled:
                self.pgrid.verifySettings()
                self.xigrid.verifySettings()
        elif self.type == TYPE_PPARPPERP:
            raise DREAMException("{}: No support implemented yet for 'ppar/pperp' grids.".format(self.name))
        else:
            raise DREAMException("{}: Unrecognized momentum grid type specified: {}.".format(self.name, self.type))


