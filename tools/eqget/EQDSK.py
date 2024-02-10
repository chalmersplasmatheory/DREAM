#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import re

from EqBase import EqBase


class EQDSK(EqBase):
    

    def __init__(self, eqdskin, cocos=1, process=True, override_psilim=False):
        """
        Constructor.

        :param eqdskin:          Name of EQDSK file to load or EQDSK structure.
        :param cocos:            COCOS number indicating sign convention of equilibrium.
        :param process:          If ``True`` calculates derived quantities from input data.
        :param override_psilim:  If ``True``, calculates the actual poloidal flux at the magnetic axis.
        """
        if type(eqdskin) == str:
            self.load(eqdskin, cocos=cocos, process=process, override_psilim=override_psilim)
        elif type(eqdskin) == dict:
            self.eqdsk = eqdskin
            if process:
                self.process_data(eqdskin, cocos=cocos, override_psilim=override_psilim)
        else:
            raise ValueError("Unrecognized input argument.")


    def load(self, filename, cocos=1, process=True, override_psilim=False):
        """
        Load and initialize this object from the named EQDSK file.
        """
        self.eqdsk = self.load_eqdsk(filename, cocos=cocos)
        if process:
            self.process_data(self.eqdsk, override_psilim=override_psilim)


    def load_eqdsk(self, filename, cocos=1):
        """
        Load EQDSK data from the specified file.
        """
        REPATTERN = r"[ +\-]?\d+(?:\.\d+(?:[Ee][\+\-]\d\d)?)?"
        eqdsk = {
            'filename': filename,
            'cocos': cocos
        }

        def _next_value(f):
            pattern = re.compile(REPATTERN)

            for line in f:
                matches = pattern.findall(line)
                for m in matches:
                    if "." in m:
                        yield float(m)
                    else:
                        yield int(m)


        def _read_1d(v, n):
            val = np.zeros(n)
            for i in range(n):
                val[i] = next(v)

            return val


        def _read_2d(v, n, m):
            val = np.zeros((n, m))
            for j in range(m):
                for i in range(n):
                    val[i,j] = next(values)

            return val


        with open(filename, 'r') as f:
            # Load header
            l = f.readline().split()
            eqdsk['stitle'] = ' '.join(l[:-3])
            eqdsk['ind1'] = int(l[-3])
            eqdsk['nx'] = int(l[-2])
            eqdsk['ny'] = int(l[-1])

            nr = eqdsk['nx']
            nz = eqdsk['ny']

            # Read meta data
            meta = []
            pattern = re.compile(REPATTERN)
            for i in range(4):
                l = f.readline()
                matches = pattern.findall(l)
                mlist = []
                for m in matches:
                    mlist.append(float(m))

                meta.append(mlist)

            values = _next_value(f)
            r1d = lambda n : _read_1d(values, n)
            r2d = lambda n, m : _read_2d(values, n, m)

            eqdsk['fpol'] = r1d(nr)
            eqdsk['p'] = r1d(nr)
            eqdsk['ffprime'] = r1d(nr)
            eqdsk['pprime'] = r1d(nr)
            eqdsk['psi'] = r2d(nr, nz)
            eqdsk['q'] = r1d(nr)

            eqdsk['psimesh'] = np.linspace(0, 1, nr)
            eqdsk['rhopsi'] = np.sqrt(eqdsk['psimesh'])

            # Load plasma boundary and limiter
            #leftover = []
            eqdsk['nbbound'] = next(values)
            eqdsk['nblim'] = next(values)
            if eqdsk['nbbound'] == int(eqdsk['nbbound']):
                # Load plasma boundary
                rzplasma = r2d(2, eqdsk['nbbound'])
                eqdsk['rplas'] = rzplasma[0,:]
                eqdsk['zplas'] = rzplasma[1,:]
            else:
                eqdsk['rplas'] = np.array([])
                eqdsk['zplas'] = np.array([])
                #leftover.append(eqdsk['nbbound'])

            if type(eqdsk['nblim']) == int:
                # Load limiter
                rzlimiter = r2d(2, eqdsk['nblim'])
                eqdsk['rlim'] = rzlimiter[0,:]
                eqdsk['zlim'] = rzlimiter[1,:]
            else:
                eqdsk['rlim'] = np.array([])
                eqdsk['zlim'] = np.array([])
                #leftover.append(eqdsk['nblim'])

            # Ignore remainder of file

        # Sort out meta data
        eqdsk['rboxlen'] = meta[0][0]
        eqdsk['zboxlen'] = meta[0][1]
        eqdsk['rcentr'] = meta[0][2]

        if eqdsk['rcentr'] < 0:
            print(f"WARNING: R0 is negative. R0 = {eqdsk['rcentr']}")

        eqdsk['rleft'] = meta[0][3]
        eqdsk['zmid'] = meta[0][4]
        eqdsk['raxis'] = meta[1][0]
        eqdsk['zaxis'] = meta[1][1]
        eqdsk['psiaxis'] = meta[1][2]
        eqdsk['psiedge'] = meta[1][3]
        eqdsk['bcentr'] = meta[1][4]
        eqdsk['ip'] = meta[2][0]

        eqdsk['R'] = np.linspace(eqdsk['rleft'], eqdsk['rleft']+eqdsk['rboxlen'], eqdsk['nx'])
        eqdsk['Z'] = np.linspace(eqdsk['zmid']-eqdsk['zboxlen']/2, eqdsk['zmid']+eqdsk['zboxlen']/2, eqdsk['ny'])

        # Construct limiter from box if needed
        if len(eqdsk['rlim']) == 1 and eqdsk['nblim'] == 1:
            eqdsk['rlim'] = [
                eqdsk['rleft'],
                eqdsk['rleft']+eqdsk['rboxlen'],
                eqdsk['rleft']+eqdsk['rboxlen'],
                eqdsk['rleft'],
                eqdsk['rleft']
            ]
        elif len(eqdsk['rlim']) == 1:
            print(f"Inconsistency: len(eqdsk['rlim']) = {len(eqdsk['rlim'])}, but eqdsk['nblim'] = {eqdsk['nblim']}.")

        return eqdsk
    

