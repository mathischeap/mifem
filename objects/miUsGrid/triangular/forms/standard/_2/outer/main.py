# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:14 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard._2.base.main import miUsTriangular_S2F_Base


class miUsTriangular_S2F_Outer(miUsTriangular_S2F_Base):
    """"""

    def __init__(self, mesh, space, name='Tri-os2f'):
        """"""
        super(miUsTriangular_S2F_Outer, self).__init__(mesh, space, 'outer', name)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_2/outer/main.py
    import numpy as np
    from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
    from tests.objects.miUsGrid.triangular.randObj.rand_mesh import mesh
    from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

    def func(t, x, y):
        return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) + t

    p = 2

    space = miUsGrid_TriangularFunctionSpace(p)

    f2 = miUsTriangular_S2F_Outer(mesh, space)

    scalar = miUsGrid_Triangular_Scalar(mesh, func)

    f2.CF = scalar

    scalar.current_time = 0

    f2.discretize()

    mesh.visualize()