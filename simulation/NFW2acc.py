import numpy as np

class NFW:
    def __init__(self, Rs, cte):
        """In code units"""
        self.Rs = Rs
        self.cte = cte

    def acc(self, pos):
        """Return acceleration due to NFW potential"""
        r = np.linalg.norm(pos, axis=1) #sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
        dPhi = -self.cte * ( np.log(1+r/self.Rs) - r/(self.Rs+r) ) / (r*r);
        return dPhi[:, np.newaxis] * pos/r[:, np.newaxis]

_Rs = 117.47806
_NFWCte = 3.2473662e+08
nfw = NFW(_Rs, _NFWCte)

if __name__ == '__main__':
    # From simulation, NFW of M=1e14 Msol (code units)
    Rs = 117.47806
    NFWCte = 3.2473662e+08
    nfw = NFW(Rs, NFWCte)