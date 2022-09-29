# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from scipy.sparse import csc_matrix
import numpy as np


class miUsGrid_TriangularFunctionSpace_IncidenceMatrix(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()

    @property
    def miUsTriangular_S0F_Outer(self):
        """curl"""
        sn = self._space_.local_numbering.miUsTriangular_S0F_Outer
        dn = self._space_.local_numbering.miUsTriangular_S1F_Outer
        E = np.zeros((self._space_.num_basis.miUsTriangular_S1F_Outer,
                      self._space_.num_basis.miUsTriangular_S0F_Outer), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                row = dn[0][i,j]
                if row >= 0:
                    E[row, sn[0][i,j]]   = -1   # L ##########
                    E[row, sn[0][i,j+1]] = +1   # R ##########
        I, J = np.shape(dn[1])
        for j in range(J):
            for i in range(I):
                row = dn[1][i,j]
                if row >= 0:
                    E[row, sn[0][i,j]]   = +1    # U
                    E[row, sn[0][i+1,j]] = -1    # D
        return csc_matrix(E)


    @property
    def miUsTriangular_S0F_Inner(self):
        """Grad"""
        sn = self._space_.local_numbering.miUsTriangular_S0F_Inner
        dn = self._space_.local_numbering.miUsTriangular_S1F_Inner
        E = np.zeros((self._space_.num_basis.miUsTriangular_S1F_Inner,
                      self._space_.num_basis.miUsTriangular_S0F_Inner), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                row = dn[0][i,j]
                if row >= 0:
                    E[row, sn[0][i,j]]   = -1   # U
                    E[row, sn[0][i+1,j]] = +1   # D
        I, J = np.shape(dn[1])
        for j in range(J):
            for i in range(I):
                row = dn[1][i,j]
                if row >= 0:
                    E[row, sn[0][i,j]]   = -1    # L
                    E[row, sn[0][i,j+1]] = +1    # R
        return csc_matrix(E)

    @property
    def miUsTriangular_S1F_Outer(self):
        """div."""
        sn = self._space_.local_numbering.miUsTriangular_S1F_Outer
        dn = self._space_.local_numbering.miUsTriangular_S2F_Outer
        E = np.zeros((self._space_.num_basis.miUsTriangular_S2F_Outer,
                      self._space_.num_basis.miUsTriangular_S1F_Outer), dtype=int)

        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                Left = sn[1][i,j]
                if Left  >= 0:
                    E[dn[0][i,j], Left] = -1    # L
                E[dn[0][i,j], sn[1][i,j+1]] = +1    # R
                E[dn[0][i,j], sn[0][i  ,j]] = -1    # U
                E[dn[0][i,j], sn[0][i+1,j]] = +1    # D
        return csc_matrix(E)

    @property
    def miUsTriangular_S1F_Inner(self):
        """row"""
        sn = self._space_.local_numbering.miUsTriangular_S1F_Inner
        dn = self._space_.local_numbering.miUsTriangular_S2F_Inner
        E = np.zeros((self._space_.num_basis.miUsTriangular_S2F_Inner,
                      self._space_.num_basis.miUsTriangular_S1F_Inner), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                Left = sn[0][i,j]
                if Left >= 0:
                    E[dn[0][i,j], Left] = +1    # L
                E[dn[0][i,j], sn[0][i,j+1]] = -1    # R
                E[dn[0][i,j], sn[1][i  ,j]] = -1    # U
                E[dn[0][i,j], sn[1][i+1,j]] = +1    # D
        return csc_matrix(E)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
