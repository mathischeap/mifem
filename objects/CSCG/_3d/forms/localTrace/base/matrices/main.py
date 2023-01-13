# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 8:29 PM
"""
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
            MM,
        )
