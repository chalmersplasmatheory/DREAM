import matplotlib.pyplot as plt
import numpy as np
from .OtherFluidQuantity import OtherFluidQuantity
from .OutputException import OutputException
from .CurrentDensity import CurrentDensity
from .ScalarQuantity import ScalarQuantity


class LCFSLoss(OtherFluidQuantity):

    def __init__(self, name, data, description, grid, output, momentumgrid):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, description=description, grid=grid, output=output)
        
    def calcScrapedOffJre(self, r=None, t=None):
    
        if r is None:
            r=slice(None)
            
        if t is None:
            t=slice(None)
            
        jreScrapedOff = self.output.eqsys.j_re.data[t,r]
        
        for i in range(len(self.grid.r[r])):
            if np.any(self.data[:,i]<0):
                itScrapeOff = np.argwhere(self.data[:,i]<0)[0][0]
                
                jreScrapedOff[itScrapeOff:,i] = jreScrapedOff[itScrapeOff,i]
        
                jreScrapedOff[:itScrapeOff,i] = 0
            else:
                jreScrapedOff[:,i] = 0
        
        return jreScrapedOff
        
    def calcScrapedOffIre(self, r=None, t=None):
    
        if r is None:
            r=slice(None)
            
        if t is None:
            t=slice(None)
            
        jreScrapedOff = CurrentDensity(name = self.name, data = self.calcScrapedOffJre(r,t), grid = self.grid, output=self.output)
        
        return jreScrapedOff.current(t = t)
        
    def plotScrapedOffJre(self, r=None, t=None, **kwargs):
        jreScrapedOff = CurrentDensity(name = 'scraped-off $J_{re}$ [A/m$^2$]', data = self.calcScrapedOffJre(), grid = self.grid, output=self.output)
        return jreScrapedOff.plot(r=r, t=t, **kwargs)
        
    def plotScrapedOffIre(self, r=None, t=None, **kwargs):
        jreScrapedOff = ScalarQuantity(name = 'scraped-off $I_{re}$ [A]', data = self.calcScrapedOffIre(), grid = self.grid, output=self.output)
        return jreScrapedOff.plot(t=t, **kwargs)
