# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 2:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from importlib import import_module


class _2nCSCG_RF2_Mesh_p_Refinement_Allocator(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID):
        """"""
        path = self.___pR_path___()[ID]
        name = self.___pR_name___()[ID]
        return getattr(import_module(path), name)(self._mesh_)

    @classmethod
    def ___pR_name___(cls):
        """"""
        return {'simple': "pR_Simple",
                }

    @classmethod
    def ___pR_path___(cls):
        """"""
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'simple': base_path + "simple",
                }


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
