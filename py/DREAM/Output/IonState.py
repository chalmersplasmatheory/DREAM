# Object representing a single ion charge state


import numpy as np
from .FluidQuantity import FluidQuantity


class IonState(FluidQuantity):
    

    def __init__(self, name, Z, Z0, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=("{}-{}".format(name, Z0)), data=data, attr=attr, grid=grid, output=output)

        self.Z  = Z
        self.Z0 = Z0


    def getRomanName(self):
        """
        Returns the name of this ion charge state using roman
        numerals for the charge state number.
        """
        val = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
        syb = ["M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"]
        roman_num = ''
        i = 0
        num = self.Z0

        if num == 0:
            roman_num = '0'
        else:
            while num > 0:
                for _ in range(num // val[i]):
                    roman_num += syb[i]
                    num -= val[i]
                i += 1

        return self.name.split('-')[0]+f'-{roman_num}'


