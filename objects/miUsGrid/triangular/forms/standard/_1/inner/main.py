# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 6:21 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard._1.base.main import miUsTriangular_S1F_Base

from objects.miUsGrid.triangular.forms.standard._1.inner.discretize.main import miUsTriangular_iS1F_Discretize
from objects.miUsGrid.triangular.forms.standard._1.inner.reconstruct.main import miUsTriangular_iS1F_Reconstruct
from objects.miUsGrid.triangular.forms.standard._1.inner.do.main import miUs_Triangular_iS1F_Do
from objects.miUsGrid.triangular.forms.standard._1.inner.matrices.main import miUs_Triangular_iS1F_Matrices
from objects.miUsGrid.triangular.forms.standard._1.inner.operators.main import miUs_Triangular_iS1F_Operators
from objects.miUsGrid.triangular.forms.standard._1.inner.IDT import miUs_Triangular_iS1F_InterfaceDofTopology
from objects.miUsGrid.triangular.forms.standard._1.inner.export.main import miUsTriangular_iS1F_Export



class miUsTriangular_S1F_Inner(miUsTriangular_S1F_Base):
    """"""

    def __init__(self, mesh, space, hybrid=False, name='Tri-is1f'):
        """"""
        super(miUsTriangular_S1F_Inner, self).__init__(mesh, space, hybrid, 'inner', name)

        self._discretize_ = miUsTriangular_iS1F_Discretize(self)
        self._reconstruct_ = miUsTriangular_iS1F_Reconstruct(self)
        self._do_ = miUs_Triangular_iS1F_Do(self)
        self._matrices_ = miUs_Triangular_iS1F_Matrices(self)
        self._operators_ = miUs_Triangular_iS1F_Operators(self)
        self._IDT_ = miUs_Triangular_iS1F_InterfaceDofTopology(self)
        self._export_ = miUsTriangular_iS1F_Export(self)

        self._freeze_self_()


    def ___Pr_check_BC_CF___(self, func):
        """"""
        assert func.__class__.__name__ in (
            "miUsGrid_Triangular_Vector", # it is considered as w, and we will w \times n on boundary
            "miUsGrid_Triangular_Scalar", # itself is considered as w \times n.
        ), f"I need a miUsGrid_Triangular_Vector or _scalar as BC.CF."
        assert func.mesh == self.mesh, f"meshes do not match!"



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/inner/main.py
    import numpy as np
    from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
    from tests.objects.miUsGrid.triangular.randObj.rand_mesh import mesh
    from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

    p = 15

    space = miUsGrid_TriangularFunctionSpace(p)

    f1 = miUsTriangular_S1F_Inner(mesh, space)

    def fx(t, x, y): return np.sin(2*np.pi*x) * np.cos(2*np.pi*y) + t
    def fy(t, x, y): return np.cos(2*np.pi*x) * np.sin(2*np.pi*y) + t

    vector = miUsGrid_Triangular_Vector(mesh, (fx, fy))

    f1.CF = vector

    vector.current_time = 0

    f1.discretize()


    print(f1.error.L())