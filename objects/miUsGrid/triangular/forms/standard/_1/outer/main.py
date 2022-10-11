# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 6:21 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard._1.base.main import miUsTriangular_S1F_Base

from objects.miUsGrid.triangular.forms.standard._1.outer.discretize.main import miUsTriangular_oS1F_Discretize
from objects.miUsGrid.triangular.forms.standard._1.outer.reconstruct.main import miUsTriangular_oS1F_Reconstruct
from objects.miUsGrid.triangular.forms.standard._1.outer.do.main import miUs_Triangular_oS1F_Do
from objects.miUsGrid.triangular.forms.standard._1.outer.matrices.main import miUs_Triangular_oS1F_Matrices
from objects.miUsGrid.triangular.forms.standard._1.outer.operators.main import miUs_Triangular_oS1F_Operators
from objects.miUsGrid.triangular.forms.standard._1.outer.IDT import miUs_Triangular_oS1F_InterfaceDofTopology
from objects.miUsGrid.triangular.forms.standard._1.outer.export.main import miUsTriangular_oS1F_Export


class miUsTriangular_S1F_Outer(miUsTriangular_S1F_Base):
    """"""

    def __init__(self, mesh, space, name='Tri-os1f'):
        """"""
        super(miUsTriangular_S1F_Outer, self).__init__(mesh, space, 'outer', name)

        self._discretize_ = miUsTriangular_oS1F_Discretize(self)
        self._reconstruct_ = miUsTriangular_oS1F_Reconstruct(self)
        self._do_ = miUs_Triangular_oS1F_Do(self)
        self._matrices_ = miUs_Triangular_oS1F_Matrices(self)
        self._operators_ = miUs_Triangular_oS1F_Operators(self)
        self._IDT_ = miUs_Triangular_oS1F_InterfaceDofTopology(self)
        self._export_ = miUsTriangular_oS1F_Export(self)

        self._freeze_self_()


    def ___Pr_check_BC_CF___(self, func):
        """"""
        assert func.__class__.__name__ in (
            "miUsGrid_Triangular_Vector", # it is considered as u, and we will use u cdot n on the boundary
            "miUsGrid_Triangular_Scalar", # itself is considered as (u cdot n) on the boundary
        ), f"I need a miUsGrid_Triangular_Vector or _scalar as BC.CF."
        assert func.mesh == self.mesh, f"meshes do not match!"

if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/outer/main.py
    import numpy as np
    from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
    from objects.miUsGrid.triangular.__test__.Random.test_mesh import mesh
    from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

    p = 15

    space = miUsGrid_TriangularFunctionSpace(p)

    f1 = miUsTriangular_S1F_Outer(mesh, space)

    def fx(t, x, y): return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t
    def fy(t, x, y): return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t


    vector = miUsGrid_Triangular_Vector(mesh, (fx, fy))

    f1.CF = vector

    vector.current_time = 0

    f1.discretize()

    print(f1.mesh)
