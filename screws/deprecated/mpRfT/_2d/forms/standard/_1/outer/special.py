# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/21 5:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix
from scipy.sparse import csr_matrix

import numpy as np

class mpRfT2_So1F_Special(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def ___Pr_test1_nonconforming_connection_through___(self, t):
        """"""
        assert t.__class__.__name__ == 'mpRfT2_NSgF'
        assert t.mesh is self._f_.mesh
        f = self._f_
        mesh = self._f_.mesh

        coo_map = mesh.coo_map.partial_segment_integral(t)

        RM = f.reconstruction.___Pr_test_matrix___.__Initialize__(coo_map)

        RM = EWC_SparseMatrix(mesh.rcfc, RM, cache_key_generator='no_cache')
        # for rc_rp in RM:
        #     if rc_rp == '1':
        #         rm = RM[rc_rp]
        #         print(rm.todense())

        return RM






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/special.py
    pass
