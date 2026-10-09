# Object representing a single ion charge state

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

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, Z=None, Z0=None):
        """
        Creates a new object of the same type where the provided quantities replace
        those of self.
        """
        if name is None:
            name = self.name
        if data is None:
            data = self.data
        if grid is None:
            grid = self.grid
        if output is None:
            output = self.output
        if attr is None:
            attr = {}
            if hasattr(self, "description"):
                attr["description"] = self.description
            if hasattr(self, "description_eqn"):
                attr["equation"] = self.description_eqn
        if Z0 is None:
            Z0 = self.Z0
        if Z is None:
            Z = self.Z
        
        return type(self)(name=name, data=data, grid=grid, output=output, attr=attr, Z=Z, Z0=Z0)
