

import numpy as np
import matplotlib.pyplot as plt

from . OutputException import OutputException


class UnknownQuantity:
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.

        name:   Name of unknown.
        data:   Data of unknown.
        attr:   List of attributes of this unknown.
        grid:   Grid used for the DREAM simulation.
        output: Parent DREAMOutput object.
        """
        self.name = name
        self.data = data
        self.grid = grid
        self.output = output

        if 'description' in attr:
            self.description = attr['description']
        if 'equation' in attr:
            self.description_eqn = attr['equation']


    def __getitem__(self, key):
        """
        Direct access to 'data' dict.
        """
        return self.data[key]


    def __add__(self, other): return self.add(other, type(self))
    def __mul__(self, other): return self.mul(other, type(self))
    def __sub__(self, other): return self.sub(other, type(self))
    def __truediv__(self, other): return self.div(other, type(self))
    def __pow__(self, other): return self.pow(other, type(self))

    def __radd__(self, other): return self.add(other, type(self), reverse=True)
    def __rmul__(self, other): return self.mul(other, type(self), reverse=True)
    def __rsub__(self, other): return self.sub(other, type(self), reverse=True)
    def __rtruediv__(self, other): return self.div(other, type(self), reverse=True)
    def __rpow__(self, other): return self.pow(other, type(self), reverse=True)


    def add(self, other, qty=None, reverse=False):
        """Add value or other UnknownQuantity to this object."""
        return self._operate(lambda a,b: a+b, other, qty, reverse, '+')


    def div(self, other, qty=None, reverse=False):
        """Divide UnknownQuantity by value or other UnknownQuantity."""
        return self._operate(lambda a,b: a/b, other, qty, reverse, '/')

    
    def mul(self, other, qty=None, reverse=False):
        """Multiply value or other UnknownQuantity to this object."""
        return self._operate(lambda a,b: a*b, other, qty, reverse, '*')


    def pow(self, other, qty=None, reverse=False):
        """Exponentiate UnknownQuantity to the given value or other UnknownQuantity."""
        return self._operate(lambda a,b: a**b, other, qty, reverse, '^')


    def sub(self, other, qty=None, reverse=False):
        """Subtract value or other UnknownQuantity to this obejct."""
        return self._operate(lambda a,b: a-b, other, qty, reverse, '-')


    def _operate(self, operation, other, qty=None, reverse=False, opname=None):
        """
        Execute the 'operation' on 'self' and 'other'. Return
        the result as a new object of type 'qty'.

        operation: Function executing the desired operation on the two
                   operands.
        other:     Second operand of addition.
        qty:       Type of object to generate (default: UnknownQuantity).
        reverse:   If 'True', reverses the order of the operands passed to
                   'operation'.
        opname:    Name of operation.
        """
        op = lambda a,b: operation(b, a) if reverse else operation(a, b)

        # Perform operation
        v = None
        otherName = 'const'
        data = self.data[:]
        if type(other) == float or type(other) == int:
            v = op(data, other)
            otherName = str(other)
        elif type(other) == np.ndarray:
            if self.data.shape != other.shape:
                raise OutputException("Mismatching dimensions of operands: {} and {}.".format(self.data.shape, other.shape))
            
            v = op(data, other)
        elif self.data.shape == other.data.shape or (np.isscalar(data) or np.isscalar(other.data)):
            v = op(data, other.data)

            # If different types, we locate the closest
            # common ancestor and convert to that type
            if type(other) != qty:
                classes = [type(self).mro(), type(other).mro()]
                for x in classes[0]:
                    if all(x in mro for mro in classes):
                        qty = x
                        break

            otherName = other.name
        else:
            raise OutputException("Unsupported type of operand: {}".format(type(other)))

        # Determine new name
        def nametransform(n):
            if ' ' in n: return '('+n+')'
            else: return n

        name1 = nametransform(self.name)
        name2 = nametransform(otherName)

        if reverse:
            name1, name2 = name2, name1

        newname = self.name
        if opname is not None:
            if type(opname) == str:
                newname = '{} {} {}'.format(name1, opname, name2)
            else:
                newname = opname(self.name, otherName)

        # Construct new object
        if qty is None:
            return UnknownQuantity(name=newname, data=v, grid=self.grid, output=self.output, attr={'description': '', 'equation': newname})
        else:
            return qty(name=newname, data=v, grid=self.grid, output=self.output, attr={'description': '', 'equation': newname})


    def getName(self): return self.name


    def getData(self): return self.data[:]


    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of ion species and
        charge states) covered by this quantity. The total number of elements
        in 'self.data' is the size of the grid on which this quantity lives
        (i.e. scalar grid, fluid grid, or a kinetic grid) times this number.
        """
        # Must be implemented specifically in derived classes
        raise NotImplementedError()
    

    def getTeXName(self):
        return self.name.replace('_', r'\_')


    def getTeXIntegralName(self):
        return 'Integrated '+self.getTeXName()


