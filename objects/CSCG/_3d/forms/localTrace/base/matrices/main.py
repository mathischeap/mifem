# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 8:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix

class _3dCSCG_LocalTrace_Matrices(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._T_ = None
        self._S_ = None
        self._freeze_self_()

    @property
    def mass(self):
        MM = self._ltf_.___PrLT_mass_matrices___()
        return EWC_SparseMatrix(
            self._ltf_.mesh,
            MM
        )

    @property
    def trace(self):
        """Return the trace matrix."""
        if self._T_ is None:
            k = self._ltf_.k
            formName = f'_3dCSCG_{int(k)}LocalTrace'
            T = getattr(self._ltf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._ltf_.mesh.elements, T, 'constant')
        return self._T_

    @property
    def selective(self):
        """Return the selective (mesh-element -> trace element) matrix.

        Like the trace matrix but without minus sign.
        """
        if self._S_ is None:
            k = self._ltf_.k
            formName = f'_3dCSCG_{int(k)}LocalTrace'
            S = getattr(self._ltf_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._ltf_.mesh.elements, S, 'constant')
        return self._S_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
