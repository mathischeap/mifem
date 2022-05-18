# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 9:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from importlib import import_module


class _2nCSCG_RF2_MeshSpacePolynomialsBasisAllocator(FrozenOnly):
    """"""

    def __init__(self, form_name):
        """"""
        NAME = self.___basis_name___()[form_name]
        PATH = self.___basis_path___()[form_name]
        self._CLASS_ = getattr(import_module(PATH), NAME)
        self._freeze_self_()

    def __call__(self, mesh, xi_eta):
        """"""
        return self._CLASS_(mesh, xi_eta)

    @classmethod
    def ___basis_name___(cls):
        return {'_2nCSCG_RF2_InnerStandard0Form': "_2nCSCG_RF2_MeshSpacePolynomialsBasis_S0F",
                '_2nCSCG_RF2_OuterStandard0Form': "_2nCSCG_RF2_MeshSpacePolynomialsBasis_S0F",
                }

    @classmethod
    def ___basis_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'_2nCSCG_RF2_InnerStandard0Form': base_path + "s0f",
                '_2nCSCG_RF2_OuterStandard0Form': base_path + "s0f",
                }



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
