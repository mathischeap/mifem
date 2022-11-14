# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 8:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

from objects.miUsGrid.triangular.fields.vector.export.vtk.main import miUsGrid_Triangular_Vector_Export_VTK


class miUsGrid_Triangular_Vector_Export(FrozenOnly):
    """"""

    def __init__(self, vf):
        """"""
        self._vf_ = vf
        self._vtk_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.vtk.triangle_and_quad(*args, **kwargs)

    @property
    def vtk(self):
        if self._vtk_ is None:
            self._vtk_ = miUsGrid_Triangular_Vector_Export_VTK(self._vf_)
        return self._vtk_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
