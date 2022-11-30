# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 10:20 AM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.mpRfT._2d.forms.standard._1.base.matrices.helpers.inner_mass import InnerMassMatrixGenerator
from objects.mpRfT._2d.forms.standard._1.base.matrices.helpers.outer_mass import OuterMassMatrixGenerator


class mpRfT2_S1F_Matrices(FrozenOnly):
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
            E = cell.space.incidence_matrix._2dCSCG_1Form_Inner
        else:
            E = cell.space.incidence_matrix._2dCSCG_1Form_Outer
        return E

    @property
    def mass(self):
        if self._f_.orientation == 'inner':
            return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                    InnerMassMatrixGenerator(self._f_),
                                    cache_key_generator=self._f_.mesh.___Pr_metric_N_key___
                                    )
        else:
            return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                    OuterMassMatrixGenerator(self._f_),
                                    cache_key_generator=self._f_.mesh.___Pr_metric_N_key___
                                    )






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/base/matrices/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    f1 = fc('1-f-o')

    E = f1.matrices.incidence

    M = f1.matrices.mass

    for rp in M:
        print(M[rp].shape)


