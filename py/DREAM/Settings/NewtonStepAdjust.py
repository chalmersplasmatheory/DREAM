# Handler for Newton step adjuster settings


TYPE_PHYSICAL = 1
TYPE_BACKTRACK = 2


class NewtonStepAdjust:
    

    def __init__(self, ttype=TYPE_PHYSICAL):
        """
        Constructor.
        """
        self.type = ttype
        self.backtrack_quantities = []


    def setBacktrack(self, *unknowns):
        """
        Apply the backtrack algorithm to the named unknown quantities. If
        no unknown quantities are named, the algorithm is applied to all
        unknowns in the equation system.
        """
        self.type = TYPE_BACKTRACK
        self.backtrack_quantities = unknowns


    def setPhysicalStep(self):
        """
        Use the maximal physical step length adjuster, which only limits
        the Newton step such that the resulting solution is physical. This
        is the default option and is applied to all other step adjustment
        algorithms.
        """
        self.type = TYPE_PHYSICAL


    def fromdict(self, data):
        """
        Load settings for this stepadjuster from the given dictionary.
        """
        self.type = int(data['type'])

        if self.type == TYPE_BACKTRACK:
            self.backtrack_quantities = data['quantities']


    def todict(self):
        """
        Convert these stepadjuster settings to a dictionary.
        """
        data = { 'type': self.type }

        if self.type == TYPE_BACKTRACK:
            data['quantities'] = ';'.join(self.backtrack_quantities)

        return data


