# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from objects.miUsGrid.triangular.forms.standard._2.base.export.vtk.main import miUsTriangular_S2F_Export_VTK


class miUsTriangular_S2F_Export(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._vtk_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.vtk.triangle_and_quad(*args, **kwargs)

    @property
    def vtk(self):
        if self._vtk_ is None:
            self._vtk_ = miUsTriangular_S2F_Export_VTK(self._sf_)
        return self._vtk_



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_2/base/export/main.py
    from __init__ import miTri
    fc = miTri.form('st32', 3)
    import numpy as np

    def p_func(t, x, y): return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t

    scalar = fc('scalar', p_func)

    f0 = fc('2-f-o', name='s2f')
    f0.CF = scalar
    scalar.current_time = 0
    f0.discretize()

    f0.export()