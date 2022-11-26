# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 9:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from components.miscellaneous.miprint import miprint
from components.miscellaneous.mirand import randint
import numpy as np
from __init__ import miTri

def _u_fun(t, x, y): return np.pi * np.exp(np.pi * x) * np.sin(np.pi * y) + 0 * t
def _v_fun(t, x, y): return np.pi * np.sin(np.pi * x) * np.cos(0.983*np.pi * y) + 0 * t



class miUsGrid_Triangles_ReconstructionMatrices(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("TriTest0 [miUsGrid_Triangles_ReconstructionMatrices] ...... ", flush=True)
        p = randint(2, 3)
        self.fc = miTri.call('rand0', p)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        f0 = self.fc('0-f-i')
        f1i = self.fc('1-f-i')
        f1o = self.fc('1-f-o')
        f2 = self.fc('2-f-o')

        scalar = self.fc('scalar', _u_fun)
        vector = self.fc('vector', [_u_fun,_v_fun])

        scalar.current_time = 0
        vector.current_time = 0

        f0.CF = scalar
        f2.CF = scalar
        f0.discretize()
        f2.discretize()
        f1i.CF = vector
        f1o.CF = vector
        f1i.discretize()
        f1o.discretize()

        xi = np.linspace(-0.99, 0.99, 7)
        et = np.linspace(-0.98, 0.98, 5)

        R0 = f0.reconstruct(xi, et, ravel=True, value_only=True)
        R2 = f2.reconstruct(xi, et, ravel=True, value_only=True)
        R1i = f1i.reconstruct(xi, et, ravel=True, value_only=True)
        R1o = f1o.reconstruct(xi, et, ravel=True, value_only=True)

        rm0 = f0.do.make_reconstruction_matrix_on_grid(xi, et)
        rm2 = f2.do.make_reconstruction_matrix_on_grid(xi, et)
        rm1i = f1i.do.make_reconstruction_matrix_on_grid(xi, et)
        rm1o = f1o.do.make_reconstruction_matrix_on_grid(xi, et)

        for e in self.fc.mesh.elements:
            r0 = rm0[e] @ f0.cochain.local[e]
            np.testing.assert_array_almost_equal(r0, R0[e][0])
            r2 = rm2[e] @ f2.cochain.local[e]
            np.testing.assert_array_almost_equal(r2, R2[e][0])
            r1i_x = rm1i[e][0] @ f1i.cochain.local[e]
            np.testing.assert_array_almost_equal(r1i_x, R1i[e][0])
            r1i_y = rm1i[e][1] @ f1i.cochain.local[e]
            np.testing.assert_array_almost_equal(r1i_y, R1i[e][1])
            r1o_x = rm1o[e][0] @ f1o.cochain.local[e]
            np.testing.assert_array_almost_equal(r1o_x, R1o[e][0])
            r1o_y = rm1o[e][1] @ f1o.cochain.local[e]
            np.testing.assert_array_almost_equal(r1o_y, R1o[e][1])

        return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/unittests/standard_forms/reconstruction_matrices.py
    miUsGrid_Triangles_ReconstructionMatrices()()
