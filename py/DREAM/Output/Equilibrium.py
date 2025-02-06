
import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import RectBivariateSpline

class Equilibrium:
    

    def __init__(self, eq=None):
        """
        Constructor.

        :param eq: Equilibrium data from DREAM output.
        """
        if eq is not None:
            self.setEquilibrium(eq)


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s  = f'TOKAMAK EQUILIBRIUM (ntheta = {self.RMinusR0.shape[0]}, npsi = {self.RMinusR0.shape[1]})\n'
        s += f'  Axis at ({self.R0[0]}, {self.Z0[0]})\n'
        s += f'  Minor radius a = {self.RMinusR0_f[0,-1]} m\n'

        return s


    def setEquilibrium(self, eq):
        """
        Set equilibrium data based on output from DREAM.
        """
        self.R0 = eq['R0']
        self.Z0 = eq['Z0']
        self.RMinusR0 = eq['RMinusR0']
        self.RMinusR0_f = eq['RMinusR0_f']
        self.ZMinusZ0 = eq['ZMinusZ0']
        self.ZMinusZ0_f = eq['ZMinusZ0_f']
        self.theta = eq['theta']


    def visualize(self, shifted=False, maxis=True, ax=None, show=None, **kwargs):
        """
        Visualize this equilibrium.

        :param shifted: If ``True``, shifts the surfaces so that (0, 0) is on the magnetic axis.
        """
        red  = (249/255, 65/255, 68/255)
        black = (87/255, 117/255, 144/255)
        gray = (190/255, 190/255, 190/255)

        # Set up axes (if not already done)

        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        rn, zn = 0, 0
        if not shifted and not np.isinf(self.R0):
            rn = self.R0
            zn = self.Z0

        # Flux surfaces
        ax.plot(self.RMinusR0+rn, self.ZMinusZ0+zn, color=gray, linewidth=1, **kwargs)
        # ...close the flux surfaces
        ax.plot(np.array([self.RMinusR0[-1,:], self.RMinusR0[0,:]])+rn, np.array([self.ZMinusZ0[-1,:], self.ZMinusZ0[0,:]])+zn, color=gray, linewidth=1, **kwargs)
        # Limiter
        ax.plot(self.RMinusR0_f[:,-1]+rn, self.ZMinusZ0_f[:,-1]+zn, color=black, linewidth=2, **kwargs)
        ax.plot(self.RMinusR0_f[(0,-1),-1]+rn, self.ZMinusZ0_f[(0,-1),-1]+zn, color=black, linewidth=2, **kwargs)
        # Magnetic axis 
        if maxis:
            ax.plot(rn, zn, 'x', color=red)
        ax.axis('equal')

        if shifted or np.isinf(self.R0):
            ax.set_xlabel('Major radius $R-R_0$ (m)')
        else:
            ax.set_xlabel('Major radius $R$ (m)')
        ax.set_ylabel('Height $Z$ (m)')

        if show:
            plt.show()

        return ax
        
        
    def getRThetaPhiFromCartesian(self, x, y, z, startingGuess = None, lengthScale=1e-2, tol=1e-4):
        """
        Calculate the flux surface coordinates corresponding to the given cartesian coordinates, centeread at the magnetic axis at phi=0 (same coordinates as used to track the SPI shards)
        
        :param x: x-coordinates for the desired points. May be an array.
        :param y: y-coordinates for the desired points. May be an array.
        :param z: z-coordinates for the desired points. May be an array.
        :param startingGuess: Starting guesses for r, which may be given to speed up the convergence. If ``None``, the toroidal minor radius coordinate rho will be used as a starting guess.
        :param lengthScale: Estimate of how far away from the starting guess the r-coordinates are likely to be. Used as a step size when looking for a range of r-values encapsulating the true values of r, to prepare for the bisection used to calculate r.
        :param tol: Determines the accuracy to which the bisection is performed. The bisection will stop when the interval becomes smaller than tol*lengthScale.
        """
        # Major radius coordinate
        R = np.hypot(x+self.R0,z)

        #Position vector
        rhox = x+self.R0 - self.R0*(x+self.R0)/R
        rhoy = y
        rhoz = z - self.R0*z/R

        # Minor radius at poloidal angle
        rho = np.sqrt(rhox**2 + rhoy**2 + rhoz**2)

        # Poloidal angle
        theta = np.zeros(R.shape)
        theta[R>=self.R0] = np.arctan2(rhoy[R>=self.R0], np.hypot(rhox[R>=self.R0],rhoz[R>=self.R0]))
        theta[R<self.R0] = np.arctan2(rhoy[R<self.R0], -np.hypot(rhox[R<self.R0], rhoz[R<self.R0]))
        theta[theta<0] = theta[theta<0]+2*np.pi


        # Bisection to find radial coordinate corresponding
        # to 'r' at 'theta'...
        # We make a guess for a valid search intervall of startingGuessR+/-lengthScale,
        # and check if it has to be expanded before actually starting with the bisection
        if startingGuess is None:
            startingGuess = rho
        ra_all = startingGuess-lengthScale
        ra_all[ra_all<0] = 0
        rb_all = ra_all + 2*lengthScale
        rho_all = rho
        theta_all = theta

        getx = RectBivariateSpline(self.RMinusR0_f[0,:], self.theta, self.RMinusR0_f[:].T)
        gety = RectBivariateSpline(self.RMinusR0_f[0,:], self.theta, self.ZMinusZ0_f[:].T)

        # Check which points are actually inside the plasma
        xxmax = getx(self.RMinusR0_f[0,-1], theta, grid=False)
        yymax = gety(self.RMinusR0_f[0,-1], theta, grid=False)
        rhomax = np.hypot(xxmax, yymax)
        inside_plasma = rho<rhomax

        ra = ra_all[inside_plasma]
        rb = rb_all[inside_plasma]
        theta = theta_all[inside_plasma]
        rho = rho_all[inside_plasma]

        # Shift ra and rb of rb is outside the grid
        ra[rb>self.RMinusR0_f[0,-1]] -= (rb[rb>self.RMinusR0_f[0,-1]] - self.RMinusR0_f[0,-1])
        ra[ra<0] = 0
        rb[rb>self.RMinusR0_f[0,-1]] -= (rb[rb>self.RMinusR0_f[0,-1]] - self.RMinusR0_f[0,-1])

        continueIntervalSearch=True

        while continueIntervalSearch: 
            xxa = getx(ra, theta, grid=False)
            yya = gety(ra, theta, grid=False)

            xxb = getx(rb, theta, grid=False)
            yyb = gety(rb, theta, grid=False)

            rhoa = np.hypot(xxa,yya)
            rhob = np.hypot(xxb,yyb)

            ra[(rhoa>rho) & (rhob>rho)] -= 2*lengthScale
            ra[ra<0] = 0
            rb = ra + 2*lengthScale

            ra[(rhoa<rho) & (rhob<rho)]+=2*lengthScale
            rb[(rhoa<rho) & (rhob<rho)]+=2*lengthScale
            
            # Shift ra and rb of rb is outside the grid
            ra[rb>self.RMinusR0_f[0,-1]] -= (rb[rb>self.RMinusR0_f[0,-1]] - self.RMinusR0_f[0,-1])
            ra[ra<0] = 0
            rb[rb>self.RMinusR0_f[0,-1]] -= (rb[rb>self.RMinusR0_f[0,-1]] - self.RMinusR0_f[0,-1])

            continueIntervalSearch = np.any(((rhoa<rho) & (rhob<rho)) | ((rhoa>rho) & (rhob>rho)))

        continueBisect=True
        while continueBisect:
            r = (ra+rb)/2
            xx = getx(r, theta, grid=False)
            yy = gety(r, theta, grid=False)
            rhor = np.hypot(xx,yy)
            ra[rhor<rho] = r[rhor<rho]
            rb[rhor>=rho] = r[rhor>=rho]
            continueBisect=np.any(np.abs(ra-rb)>lengthScale*tol)

        rb_all[inside_plasma] = r
        rb_all[~inside_plasma] = self.RMinusR0_f[0,-1]+1e-2 # Arbitrary value outside the plasma

        phi = np.arctan2(z,self.R0+x)

        return rb_all, theta, phi


