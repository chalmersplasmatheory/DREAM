import numpy as np
from scipy import interpolate

    
class geqdsk():

    def __init__(self, filename):

        def splitter(inv, step=16):
            value = []
            for k in range(int(len(inv) / step)):
                value.append(inv[step * k:step * (k + 1)])
            return value

        def merge(inv):
            return ''.join(inv)

        self.filename = filename

        with open(self.filename, 'r') as f:		EQDSK = f.read().splitlines()

        self.Ginfo = {}

        # first line is description and sizes
        self.Ginfo['CASE'] = np.array(splitter(EQDSK[0][0:48], 8))
        try:
            tmp = list([_f for _f in EQDSK[0][48:].split(' ') if _f])
            [IDUM, self.Ginfo['NW'], self.Ginfo['NH']] = list(map(int, tmp[:3]))
        except ValueError: # Can happen if no space between numbers, such as 10231023
            try:
            	IDUM = int(EQDSK[0][48:52])
            	self.Ginfo['NW'] = int(EQDSK[0][52:56])
            	self.Ginfo['NH'] = int(EQDSK[0][56:60])
            	tmp = []
            	#printd('IDUM, NW, NH',IDUM,self.Ginfo['NW'],self.Ginfo['NH'],topic='OMFITgeqdsk.load')
            except ValueError: 
            	tmp = list([_f for _f in EQDSK[0][40:].split(' ') if _f])
            	[IDUM, self.Ginfo['NW'], self.Ginfo['NH']] = list(map(int, tmp[:3]))   
        if len(tmp) > 3:
            self.Ginfo['EXTRA_HEADER'] = EQDSK[0][49 + len(re.findall('%d +%d +%d ' % (IDUM, self.Ginfo['NW'], self.Ginfo['NH']), EQDSK[0][49:])[0]) + 2:]
        offset = 1

        # now, the next 20 numbers (5 per row)
        [self.Ginfo['RDIM'], self.Ginfo['ZDIM'], self.Ginfo['RCENTR'], self.Ginfo['RLEFT'], self.Ginfo['ZMID'], \
         self.Ginfo['RMAXIS'], self.Ginfo['ZMAXIS'], self.Ginfo['SIMAG'], self.Ginfo['SIBRY'], self.Ginfo['BCENTR'], \
         self.Ginfo['CURRENT'], self.Ginfo['SIMAG'], XDUM, self.Ginfo['RMAXIS'], XDUM, \
         self.Ginfo['ZMAXIS'], XDUM, self.Ginfo['SIBRY'], XDUM, XDUM] = list(map(eval, splitter(merge(EQDSK[offset:offset + 4]))))
        offset = offset + 4

        # now I have to read NW elements
        nlNW = int(np.ceil(self.Ginfo['NW'] / 5.))
        self.Ginfo['FPOL'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self.Ginfo['PRES'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self.Ginfo['FFPRIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        self.Ginfo['PPRIME'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW
        try:
            # official gEQDSK file format saves PSIRZ as a single flat array of size rowsXcols
            nlNWNH = int(np.ceil(self.Ginfo['NW'] * self.Ginfo['NH'] / 5.))
            self.Ginfo['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset + nlNWNH])), dtype=np.float64),(self.Ginfo['NH'], self.Ginfo['NW']))
            offset = offset + nlNWNH
        except ValueError:
            # sometimes gEQDSK files save row by row of the PSIRZ grid (eg. FIESTA code)
            nlNWNH = self.Ginfo['NH'] * nlNW
            self.Ginfo['PSIRZ'] = np.reshape(np.fromiter(splitter(''.join(EQDSK[offset:offset + nlNWNH])), dtype=np.float64),(self.Ginfo['NH'], self.Ginfo['NW']))
            offset = offset + nlNWNH
        self.Ginfo['QPSI'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
        offset = offset + nlNW

        # now vacuum vessel and limiters
        [self.Ginfo['NBBBS'], self.Ginfo['LIMITR']] = list(map(int, [_f for _f in EQDSK[offset:offset + 1][0].split(' ') if _f]))
        offset = offset + 1

        nlNBBBS = int(np.ceil(self.Ginfo['NBBBS'] * 2 / 5.))
        self.Ginfo['RBBBS'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNBBBS]))))[0::2])
        self.Ginfo['ZBBBS'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNBBBS]))))[1::2])
        offset = offset + max(nlNBBBS, 1)

        try:
            # this try/except is to handle some gEQDSK files written by older versions of ONETWO
            nlLIMITR = int(np.ceil(self.Ginfo['LIMITR'] * 2 / 5.))
            self.Ginfo['RLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[0::2])
            self.Ginfo['ZLIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlLIMITR]))))[1::2])
            offset = offset + nlLIMITR
        except ValueError:
            # if it fails make the limiter as a rectangle around the plasma boundary that does not exceed the computational domain
            self.Ginfo['LIMITR'] = 5
            dd = self.Ginfo['RDIM'] / 10.
            R = np.linspace(0, self.Ginfo['RDIM'], 2) + self.Ginfo['RLEFT']
            Z = np.linspace(0, self.Ginfo['ZDIM'], 2) - self.Ginfo['ZDIM'] / 2. + self.Ginfo['ZMID']
            self.Ginfo['RLIM'] = np.array([max([R[0], np.min(self.Ginfo['RBBBS']) - dd]),
                                        min([R[1], np.max(self.Ginfo['RBBBS']) + dd]),
                                        min([R[1], np.max(self.Ginfo['RBBBS']) + dd]),
                                        max([R[0], np.min(self.Ginfo['RBBBS']) - dd]),
                                        max([R[0], np.min(self.Ginfo['RBBBS']) - dd])])
            self.Ginfo['ZLIM'] = np.array([max([Z[0], np.min(self.Ginfo['ZBBBS']) - dd]),
                                        max([Z[0], np.min(self.Ginfo['ZBBBS']) - dd]),
                                        min([Z[1], np.max(self.Ginfo['ZBBBS']) + dd]),
                                        min([Z[1], np.max(self.Ginfo['ZBBBS']) + dd]),
                                        max([Z[0], np.min(self.Ginfo['ZBBBS']) - dd])])

        try:
            [self.Ginfo['KVTOR'], self.Ginfo['RVTOR'], self.Ginfo['NMASS']] = list(map(float,[_f for _f in EQDSK[offset:offset + 1][0].split(' ') if _f]))
            offset = offset + 1

            if self.Ginfo['KVTOR'] > 0:
                self.Ginfo['PRESSW'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
                offset = offset + nlNW
                self.Ginfo['PWPRIM'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
                offset = offset + nlNW

            if self.Ginfo['NMASS'] > 0:
                self.Ginfo['DMION'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
                offset = offset + nlNW

            self.Ginfo['RHOVN'] = np.array(list(map(float, splitter(merge(EQDSK[offset:offset + nlNW])))))
            offset = offset + nlNW
        except Exception:
            pass

         # ------------------ Other quantities

        self.Ginfo['R'] = np.linspace(0,self.Ginfo['RDIM'], self.Ginfo['NW']) + self.Ginfo['RLEFT']
        self.Ginfo['Z'] = np.linspace(0,self.Ginfo['ZDIM'], self.Ginfo['NH']) - self.Ginfo['ZDIM']/2. + self.Ginfo['ZMID']

        #poloidal flux and normalized poloidal flux
        self.Ginfo['PSI']      = np.linspace(self.Ginfo['SIMAG'], self.Ginfo['SIBRY'], len(self.Ginfo['PRES']))
        self.Ginfo['PSI_NORM'] = np.linspace(0., 1., len(self.Ginfo['PRES']))

        self.Ginfo['PSIRZ']      = self.Ginfo['PSIRZ']
        self.Ginfo['PSIRZ_NORM'] = abs((self.Ginfo['PSIRZ']-self.Ginfo['SIMAG'])/(self.Ginfo['SIBRY']-self.Ginfo['SIMAG']))
        #rho poloidal
        self.Ginfo['RHOp']   = np.sqrt(self.Ginfo['PSI_NORM'])
        self.Ginfo['RHOpRZ'] = np.sqrt(self.Ginfo['PSIRZ_NORM'])

        #extend functions in PSI to be clamped at edge value when outside of PSI range (i.e. outside of LCFS)
        dp = self.Ginfo['PSI'][1]-self.Ginfo['PSI'][0]
        ext_psi_mesh = np.hstack((self.Ginfo['PSI'][0]-dp*1E6, self.Ginfo['PSI'], self.Ginfo['PSI'][-1]+dp*1E6))
        def ext_arr(inv):
            return np.hstack((inv[0], inv, inv[-1]))
        #map functions in PSI to RZ coordinate
        for name1 in ['FPOL','PRES','QPSI','FFPRIM','PPRIME']:
            self.Ginfo[name1+'RZ']=interpolate.interp1d(ext_psi_mesh, ext_arr(self.Ginfo[name1]), kind='linear', bounds_error=False)(self.Ginfo['PSIRZ'])
