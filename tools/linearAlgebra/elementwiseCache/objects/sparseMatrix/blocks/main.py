# -*- coding: utf-8 -*-
"""
If we use bmat to form EWC_SparseMatrix. We can get some blocks of the EWC_SparseMatrix.

"""


import numpy as np
from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.operators.bmat.main import bmat


class EWC_SpaMat_Blocks(FrozenOnly):
    """"""
    def __init__(self, spa_mat):
        if spa_mat.bmat_shape is False:
            raise Exception(f"Blocks only valid for bmat EWC_SparseMatrix.")
        GM0, GM1 = spa_mat.gathering_matrices
        assert GM0 is not None and GM1 is not None, f"To access blocks, must have gathering matrices."

        self._shape_ = spa_mat.bmat_shape
        self.___BLOCKS___ = np.array(spa_mat._DG_.blocks, dtype=object)


        #-- replace all None blocks by zero-blocks --------------------------------------------
        SHAPE = self.___BLOCKS___.shape
        I, J = SHAPE
        for i in range(I):
            for j in range(J):
                if self.___BLOCKS___[i, j] is None:
                    num_basis_0 = GM0.GMs[i].GLOBAL_shape[1]
                    num_basis_j = GM1.GMs[j].GLOBAL_shape[1]
                    self.___BLOCKS___[i, j] = spa_mat.__class__(
                        spa_mat.elements,(num_basis_0, num_basis_j))

        self._spa_mat_ = spa_mat
        self._freeze_self_()

    @property
    def shape(self):
        return self._shape_

    def __getitem__(self, item):
        B = self.___BLOCKS___[item]

        if B.__class__.__name__ == 'EWC_SparseMatrix':
            return B
        else:
            B = bmat(B)
            return B
