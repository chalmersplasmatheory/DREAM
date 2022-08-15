# This class provides routines for reading input profiles from CORSICA data files

import numpy as np

class ITER_Profile_reader:

    def __init__(self, filename):
        self.filename = filename
        
        f = open(self.filename,'r')
        lines = f.readlines()
        iLine = 0
        keep_reading = True
        
        # Go through the file and find basic scalar data 
        # (number of columns ncol for the profile datasets, 
        # number of radial gridpoints N_Rho and the plasma current Ip)
        while keep_reading:
            line = lines[iLine]
            #print(line)
            if line.split():
                if line.split()[0] == 'ncol':
                    self.ncol = int(lines[iLine+1])
                    iLine+=1
                if line.split()[0] == 'N_Rho':
                    self.N_Rho = int(lines[iLine+1])
                    iLine+=1
                    
                if line.split()[0] == 'plasma' and line.split()[1] == 'current':
                    self.Ip = float(lines[iLine+1].replace('D','e'))
                    iLine+=1
                            
            
            iLine+=1        
            keep_reading = iLine<len(lines)
        
    def load_profile(self, quantity):
        # This function goes through the file and finds the profile dataset named 'quantity' 
        # in the input file, where 'quantity' should be a string matching one of the 
        # quantity names in the CORSICA file, and returns this profile as a numpy array
        
        f = open(self.filename,'r')
        lines = f.readlines()
        iLine = 0
        keep_reading = True
        while keep_reading:
            line = lines[iLine]
            if line.split():
                if line.split()[0] == quantity:
                    values = np.zeros(self.N_Rho)
                    iRho = 0
                    while iRho<self.N_Rho:
                        line = lines[iLine+1]
                        for iCol in range(len(line.split())):
                            values[iRho] = float(line.split()[iCol].replace('D','e'))
                            iRho+=1
                        iLine+=1
                            
            
            iLine+=1        
            keep_reading = iLine<len(lines)
        return values
        
    def get_Ip(self):
        # Getter function for the plasma current
        return self.Ip
