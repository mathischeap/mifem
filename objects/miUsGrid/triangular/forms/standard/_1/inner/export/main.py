# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from objects.miUsGrid.triangular.forms.standard._1.inner.export.vtk.main import miUsTriangular_iS1F_Export_VTK




class miUsTriangular_iS1F_Export(FrozenOnly):
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
            self._vtk_ = miUsTriangular_iS1F_Export_VTK(self._sf_)
        return self._vtk_



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/inner/export/main.py
    from __init__ import miTri
    import numpy as np

    fc = miTri.form('st16', 2)

    def p_func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t
    def q_func(t, x, y): return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) + t

    v = fc('vector', [p_func,q_func])

    f1 = fc('1-f-i')
    f1.CF = v
    v.current_time = 0
    f1.discretize()

    f1.export()


