# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/23 5:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.miscellaneous.miprint import miprint
import numpy as np

from tests.objects.miUsGrid.triangular.randObj.test_mesh import mesh as tm
from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace
from objects.miUsGrid.triangular.forms.standard._0.inner.main import miUsTriangular_S0F_Inner
from objects.miUsGrid.triangular.forms.standard._0.outer.main import miUsTriangular_S0F_Outer
from objects.miUsGrid.triangular.forms.standard._1.inner.main import miUsTriangular_S1F_Inner
from objects.miUsGrid.triangular.forms.standard._1.outer.main import miUsTriangular_S1F_Outer

from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector

class miUsGrid_Triangle_Incidence_matrices(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("TriTest0 [miUsGrid_Triangle_Incidence_matrices] ...... ", flush=True)
        self.mesh = tm
        self.space0 = miUsGrid_TriangularFunctionSpace(1)
        self.space1 = miUsGrid_TriangularFunctionSpace(3)
        self.space2 = miUsGrid_TriangularFunctionSpace(5)
        self.space3 = miUsGrid_TriangularFunctionSpace(7)

        self.scalar = miUsGrid_Triangular_Scalar(self.mesh, self.p)
        self.vector = miUsGrid_Triangular_Vector(self.mesh, [self.fx, self.fy])

        self.scalar.current_time = 0
        self.vector.current_time = 0

        self.grad_scalar = self.scalar.numerical.grad
        self.curl_scalar = self.scalar.numerical.curl
        self.grad_scalar.current_time = 0
        self.curl_scalar.current_time = 0

        self.div_vector = self.vector.numerical.div
        self.rot_vector = self.vector.numerical.rot
        self.div_vector.current_time = 0
        self.rot_vector.current_time = 0

        self._freeze_self_()

    @staticmethod
    def p(t, x, y):
        return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + t

    @staticmethod
    def fx(t, x, y):
        return np.cos(2*np.pi*x) * np.sin(3*np.pi*y) + t

    @staticmethod
    def fy(t, x, y):
        return np.sin(2*np.pi*x) * np.cos(1.5*np.pi*y) + t

    def __call__(self):
        """"""
        Dr0 = list()
        Dr1 = list()
        Dr2 = list()
        Dr3 = list()

        for space in (self.space0, self.space1, self.space2, self.space3):
            f0i = miUsTriangular_S0F_Inner(self.mesh, space)
            f0o = miUsTriangular_S0F_Outer(self.mesh, space)
            f1i = miUsTriangular_S1F_Inner(self.mesh, space)
            f1o = miUsTriangular_S1F_Outer(self.mesh, space)

            f0i.CF = self.scalar
            f0i.discretize()

            f0o.CF = self.scalar
            f0o.discretize()

            f1i.CF = self.vector
            f1i.discretize()

            f1o.CF = self.vector
            f1o.discretize()

            df0i = f0i.coboundary()
            df0i.CF = self.grad_scalar
            dr = df0i.error.L()
            Dr0.append(dr)

            df0o = f0o.coboundary()
            df0o.CF = self.curl_scalar
            dr = df0o.error.L()
            Dr1.append(dr)

            df1i = f1i.coboundary()
            df1i.CF = self.rot_vector
            dr = df1i.error.L()
            Dr2.append(dr)

            df1o = f1o.coboundary()
            df1o.CF = self.div_vector
            dr = df1o.error.L()
            Dr3.append(dr)

        Dr0 = np.array(Dr0)
        Dr1 = np.array(Dr1)
        Dr2 = np.array(Dr2)
        Dr3 = np.array(Dr3)

        assert all(np.diff(Dr0) < 0)
        assert all(np.diff(Dr1) < 0)
        assert all(np.diff(Dr2) < 0)
        assert all(np.diff(Dr3) < 0)

        assert Dr0[-1] < 0.01
        assert Dr1[-1] < 0.01
        assert Dr2[-1] < 0.06
        assert Dr3[-1] < 0.04

        return 1


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/unittests/standard_forms/incidence_matrices.py
    miUsGrid_Triangle_Incidence_matrices()()
