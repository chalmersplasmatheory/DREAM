

from .StaticQuantity import StaticQuantity


class StaticQuantities:
    

    SPECIAL_TREATMENT = {
        # List of static quantities with their own classes
        # (as in OtherQuantities.py)
    }


    def __init__(self, name, static=None):
        """
        Constructor.
        """
        self.name = name
        self.quantities = {}
        
        if static is not None:
            self.setQuantities(static)


    def __contains__(self, index):
        return (index in self.quantities)

    
    def __getitem__(self, index):
        """
        Direct access by name to the list of quantities.
        """
        return self.quantities[index]


    def __iter__(self):
        """
        Iterate over static quantities.
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
        Get a list with the names of quantities stored in the
        output file.
        """
        return list(self.quantities.keys())


    def setQuantity(self, name, data, attributes={}):
        """
        Add the given quantity to the list of static quantities.

        name: Name of the quantity to be added.
        data: Data of the quantity to be added.
        attributes: Optional dictionary with attributes of the quantity.
        """
        desc = ""
        if 'description' in attributes:
            desc = attributes['description']

        if name in self.SPECIAL_TREATMENT:
            o = self.SPECIAL_TREATMENT[name](name=name, data=data, description=desc)
        else:
            o = StaticQuantity(name=name, data=data, description=desc)

        setattr(self, name, o)
        self.quantities[name] = o
    

    def setQuantities(self, quantities):
        """
        Add a list of static quantities to this handler.
        """
        for qty in quantities:
            # Skip attribute containers
            if qty[-2:] == '@@': continue

            # Is there an attribute container for this quantity?
            if qty+'@@' in quantities:
                self.setQuantity(name=qty, data=quantities[qty], attributes=quantities[qty+'@@'])
            else:
                self.setQuantity(name=qty, data=quantities[qty])


    def tostring(self, padding=''):
        """
        Convert this object to a string.
        """
        s = ""
        for qty in self.quantities:
            q = self.quantities[qty]

            s += f"{padding}{q.name:28s} -- {q.description}\n"

        return s


