# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 2:58 PM
"""
import sys

if './' not in sys.path: sys.path.append('../')

from screws.freeze.base import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from objects.mpRfT._2d.forms.standard._2.base.matrices.helpers.mass import MassMatrixGenerator


class mpRfT2_S2F_Matrices(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    @property
    def mass(self):
        return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                MassMatrixGenerator(self._f_),
                                cache_key_generator=self._f_.mesh.___Pr_metric_N_key___
                                )






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
