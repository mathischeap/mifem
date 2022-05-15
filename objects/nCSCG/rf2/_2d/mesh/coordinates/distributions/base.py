# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 2:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MRF2_CooDistributionBase(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._distribution_ = None

    @property
    def distribution(self):
        return self._distribution_

    @property
    def ___Pr_is_2nCSCG_RF2_mesh_coo___(self):
        return True


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
