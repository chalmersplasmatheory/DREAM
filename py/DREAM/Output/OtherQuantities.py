
import numpy as np

from . OtherQuantity import OtherQuantity
from . OtherFluidQuantity import OtherFluidQuantity
from . OtherIonSpeciesKineticQuantity import OtherIonSpeciesKineticQuantity
from . OtherKineticQuantity import OtherKineticQuantity
from . OtherScalarQuantity import OtherScalarQuantity

from . AvalancheGrowthRate import AvalancheGrowthRate


class OtherQuantities:
    

    SPECIAL_TREATMENT = {
        # List of other quantities with their own classes
        'f_hot_ripple_pmn': OtherQuantity,
        'f_re_ripple_pmn': OtherQuantity,
        'kinioniz_vsigma': OtherIonSpeciesKineticQuantity,
        'GammaAva': AvalancheGrowthRate,
        'nu_D_f1': OtherKineticQuantity,
        'nu_D_f2': OtherKineticQuantity,
        'nu_s_f1': OtherKineticQuantity,
        'nu_s_f2': OtherKineticQuantity,
    }


    def __init__(self, name, other=None, grid=None, output=None, momentumgrid=None):
        """
        Constructor.
        """
        self.name = name
        self.grid = grid
        self.quantities = {}
        self.output = output
        self.momentumgrid = momentumgrid
        
        if other is not None:
            self.setQuantities(other)


    def __contains__(self, index):
        return (index in self.quantities)


    def __getitem__(self, index):
        """
        Direct access by name to the list of quantities.
        """
        return self.quantities[index]


    def __iter__(self):
        """
        Iterate over other quantities.
        """
        for key, item in self.quantities.items():
            yield key, item


    def __repr__(self):
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return self.tostring()


    def keys(self): return self.quantities.keys()


    def getQuantityNames(self):
        """
        Get a list with the names of all other quantities
        stored in the output file.
        """
        return list(self.quantities.keys())


    def setGrid(self, grid):
        """
        Sets the grid that was used for the DREAM simulation.
        """
        self.grid = grid


    def setQuantity(self, name, data, attributes=None, datatype=None):
        """
        Add the given quantity to the list of other quantities.

        name: Name of the quantity.
        data: Data of the quantity (raw, as a dict from the output file).
        """
        desc = ""
        if 'description' in attributes:
            desc = attributes['description']

        if datatype is not None:
            if data.ndim == 4 and self.momentumgrid is not None:
                o = datatype(name=name, data=data, description=desc, grid=self.grid, output=self.output)
            else:
                o = datatype(name=name, data=data, description=desc, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid)
        elif name in self.SPECIAL_TREATMENT:
            if data.ndim == 5 and self.momentumgrid is not None:
                o = self.SPECIAL_TREATMENT[name](name=name, data=data, description=desc, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid)
            elif data.ndim == 4 and self.momentumgrid is not None:
                o = self.SPECIAL_TREATMENT[name](name=name, data=data, description=desc, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid)
            else:
                o = self.SPECIAL_TREATMENT[name](name=name, data=data, description=desc, grid=self.grid, output=self.output)
        else:
            if self.name == 'scalar':
                o = OtherScalarQuantity(name=name, data=data, description=desc, grid=self.grid, output=self.output)
            elif data.ndim == 2:
                o = OtherFluidQuantity(name=name, data=data, description=desc, grid=self.grid, output=self.output)
            elif data.ndim == 4 and self.momentumgrid is not None:
                o = OtherKineticQuantity(name=name, data=data, description=desc, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid)
            else:
                #raise Exception("Unrecognized number of dimensions of other quantity '{}': {}.".format(name, data.ndim))
                o = OtherQuantity(name=name, data=data, description=desc, grid=self.grid, output=self.output)

        setattr(self, name, o)
        self.quantities[name] = o


    def setQuantities(self, quantities):
        """
        Add a list of other quantities to this handler.
        """
        for oqn in quantities:
            # Skip attribute containers
            if oqn[-2:] == '@@': continue

            # Is there an attribute container for this quantity?
            if oqn+'@@' in quantities:
                self.setQuantity(name=oqn, data=quantities[oqn], attributes=quantities[oqn+'@@'])
            else:
                self.setQuantity(name=oqn, data=quantities[oqn])


    def resetQuantity(self, quantity, datatype):
        """
        Reinitializes the named other quantity, making it of the given
        data type.
        """
        if quantity not in self.quantities:
            return

        attr = {}
        q = self.quantities[quantity]
        if hasattr(q, 'description'):
            attr['description'] = q.description

        self.setQuantity(name=quantity, data=q.data, attributes=attr, datatype=datatype)


    def tostring(self, padding=''):
        """
        Convert this object to a string.
        """
        s = ""
        for oqn in self.quantities:
            q = self.quantities[oqn]

            s += "{}{:28s} -- {}\n".format(padding, q.name, q.description)
        
        return s


