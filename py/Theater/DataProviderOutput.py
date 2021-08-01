# Implementation of a data provider encapsulating a DREAMOutput object.


from . DataProvider import DataProvider


class DataProviderOutput(DataProvider):
    

    def __init__(self, output):
        """
        Constructor.

        :param output: DREAMOutput object to encapsulate.
        """
        self.output = output


    def getOtherInfo(self, name=None):
        """
        Returns basic information about the unknowns of
        the equation system.
        """
        if name is None:
            def addQuantities(category, oqs):
                if oqs is None:
                    return {}

                d = {}
                for oqname, oq in oqs:
                    d['{}/{}'.format(category, oqname)] = {
                        'description': oq.description,
                        'nelements': oq.data.size,
                        'nmultiples': oq.getMultiples()
                    }

                return d
                    

            other = self.output.other
            d = {}
            d = {**d, **addQuantities('fluid', other.fluid)}
            d = {**d, **addQuantities('hottail', other.hottail)}
            d = {**d, **addQuantities('runaway', other.runaway)}
            d = {**d, **addQuantities('scalar', other.scalar)}

            return d
        else:
            cat, _, name = name.partition('/')
            oq = self.output.other[cat][name]

            d = {
                'description': oq.description,
                'nelements': oq.data.size,
                'nmultiples': oq.getMultiples()
            }

            return d


    def getUnknownInfo(self, name=None):
        """
        Returns basic information about the unknowns of the
        equation system.
        """
        d = {}

        def getInfo(name):
            u = self.output.eqsys[name]
            return {
                'description': u.description,
                'equation': u.description_eqn,
                'nelements': u.data.size,
                'nmultiples': u.getMultiples()
            }

        if name is None:
            for uqn in self.output.eqsys.unknowns:
                d[uqn] = getInfo(uqn)
            return d
        else:
            return getInfo(name)


