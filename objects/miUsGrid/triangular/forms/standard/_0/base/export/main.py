# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:52 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly
from objects.miUsGrid.triangular.forms.standard._0.base.export.vtk.main import miUsTriangular_S0F_Export_VTK


class miUsTriangular_S0F_Export(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._vtk_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.vtk.triangle_and_quad(*args, **kwargs)

    @property
    def vtk(self):
        if self._vtk_ is None:
            self._vtk_ = miUsTriangular_S0F_Export_VTK(self._sf_)
        return self._vtk_


if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_0/base/export/main.py
    from __init__ import miTri
    import numpy as np

    fc = miTri.call('st16', 2)

    def p_func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t

    scalar = fc('scalar', p_func)

    f0 = fc('0-f-o')
    f0.CF = scalar
    scalar.current_time = 0
    f0.discretize()

    f0.export()
