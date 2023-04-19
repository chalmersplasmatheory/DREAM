# This class provides routines for reading input profiles from CORSICA data files

import numpy as np

class CORSICAProfileReader:

    def __init__(self, filename):
        self.filename = filename
        
        self.ncol = None
        self.N_Rho = None
        self.Ip = None
        
        self.load_basics()     
            
    def load_basics(self):  
        # This function goes through the file and finds basic scalar data 
        # (number of columns ncol for the profile datasets, 
        # number of radial gridpoints N_Rho and the plasma current Ip)
        
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        iLine = 0
        
        while iLine<len(lines):
            line = lines[iLine]
            line_content = line.split()
            if line_content:
                if line_content[0] == 'ncol':
                    self.ncol = int(lines[iLine+1])
                    iLine+=1
                if line_content[0] == 'N_Rho':
                    self.N_Rho = int(lines[iLine+1])
                    iLine+=1
                    
                if line_content[0] == 'plasma' and line_content[1] == 'current':
                    self.Ip = float(lines[iLine+1].replace('D','e'))
                    iLine+=1
                            
            
            iLine+=1     
        
    def load_profile(self, quantity):
        # This function goes through the file and finds the profile dataset named 'quantity' 
        # in the input file, where 'quantity' should be a string matching one of the 
        # quantity names in the CORSICA file, and returns this profile as a numpy array
        
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        iLine = 0
        while iLine<len(lines):
            line = lines[iLine]
            line_content = line.split()
            if line_content:
                if line_content[0] == quantity:
                    values = np.zeros(self.N_Rho)
                    iRho = 0
                    while iRho<self.N_Rho:
                        line = lines[iLine+1]
                        line_content = line.split()
                        for iCol in range(len(line_content)):
                            values[iRho] = float(line_content[iCol].replace('D','e'))
                            iRho+=1
                        iLine+=1
                            
            
            iLine+=1        

        return values
        
    def get_Ip(self):
        # Getter function for the plasma current
        return self.Ip
