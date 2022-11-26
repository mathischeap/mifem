# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import numpy as np


class miUsGrid_TriangularFunctionSpace_LocalNumbering(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()

    @property
    def miUsTriangular_S0F_Outer(self):
        p = self._space_.p
        ln = np.arange( (p+1)**2 ).reshape((p+1, p+1), order='F')
        ln -= p
        ln[:,0] = 0
        return ln, # do not remove comma

    @property
    def miUsTriangular_S0F_Inner(self):
        p = self._space_.p
        ln = np.arange( (p+1)**2 ).reshape((p+1, p+1), order='F')
        ln -= p
        ln[:,0] = 0
        return ln, # do not remove comma

    @property
    def miUsTriangular_S1F_Outer(self):
        p = self._space_.p
        LN_dy = np.arange(0, p*(p+1)).reshape((p+1, p), order='F')

        LN_dx = np.arange(p*(p+1), 2*p*(p+1)).reshape((p, p+1), order='F')

        LN_dx -= p

        LN_dx[:,0] = -1

        return LN_dy, LN_dx

    @property
    def miUsTriangular_S1F_Inner(self):
        p = self._space_.p
        LN_dx = np.arange(0, p*(p+1)).reshape((p, p+1), order='F')

        LN_dx -= p

        LN_dx[:, 0] = -1

        LN_dy = np.arange(p * p, p * p + p * (p+1)).reshape((p+1, p), order='F')

        return LN_dx, LN_dy

    @property
    def miUsTriangular_S2F_Outer(self):
        p = self._space_.p
        ln = np.arange( p**2 ).reshape((p, p), order='F')
        return ln, # do not remove comma

    @property
    def miUsTriangular_S2F_Inner(self):
        p = self._space_.p
        ln = np.arange( p**2 ).reshape((p, p), order='F')
        return ln, # do not remove comma




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
