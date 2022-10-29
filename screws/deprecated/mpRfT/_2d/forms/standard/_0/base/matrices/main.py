# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/23 9:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.mpRfT._2d.forms.standard._0.base.matrices.helpers.mass import MassMatrixGenerator

class mpRfT2_S0F_Matrices(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    @property
    def incidence(self):
        return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                self.___Pr_incidence___,
                                cache_key_generator=self._f_.mesh.___Pr_N_key___
                                )

    def ___Pr_incidence___(self, rp):
        """"""
        cell = self._f_.mesh[rp]
        if self._f_.orientation == 'inner':
            im = cell.space.incidence_matrix._2dCSCG_0Form_Inner
        else:
            im = cell.space.incidence_matrix._2dCSCG_0Form_Outer
        return im

    @property
    def mass(self):
        return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                MassMatrixGenerator(self._f_),
                                cache_key_generator=self._f_.mesh.___Pr_metric_N_key___
                                )




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/base/matrices/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    f0 = fc('0-f-o')

    E = f0.matrices.incidence