# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 3:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsGrid_TriangularMesh_Element_Visualize(FrozenOnly):
    """"""

    def __init__(self, element):
        """"""
        self._element_ = element
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""

    def matplot(self):
        """"""


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
