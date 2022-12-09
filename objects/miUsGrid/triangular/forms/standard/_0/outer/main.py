# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:14 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard._0.base.main import miUsTriangular_S0F_Base


class miUsTriangular_S0F_Outer(miUsTriangular_S0F_Base):
    """"""

    def __init__(self, mesh, space, hybrid=False, name='Tri-os0f'):
        """"""
        super(miUsTriangular_S0F_Outer, self).__init__(mesh, space, hybrid, 'outer', name)
        self._freeze_self_()






if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_0/outer/main.py
    import numpy as np

    from tests.objects.miUsGrid.triangular.randObj.test_mesh import mesh
    from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace
    from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar

    # mesh.visualize()

    p = 3

    space = miUsGrid_TriangularFunctionSpace(p)


    s0f = miUsTriangular_S0F_Outer(mesh, space)

    def func(t, x, y): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + t

    scalar = miUsGrid_Triangular_Scalar(mesh, func)

    s0f.CF = scalar

    scalar.current_time = 0

    s0f.discretize()

    # print(s0f.cochain.local)

    GM = s0f.numbering.gathering

    print(GM)

    # error = s0f.error.L()

