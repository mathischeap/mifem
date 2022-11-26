# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 11:37 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from components.miscellaneous.miprint import miprint
import numpy as np
from __init__ import miTri

def _u_fun(t, x, y): return np.pi * np.exp(np.pi * x) * np.sin(np.pi * y) + 0 * t
def _v_fun(t, x, y): return np.pi * np.sin(np.pi * x) * np.cos(0.983*np.pi * y) + 0 * t


class miUsGrid_Triangles_MassMatrices(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("TriMassMat [miUsGrid_Triangles_MassMatrices] ...... ", flush=True)
        self.fc = miTri.call('rand0', 3)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        f0 = self.fc('0-f-i')
        f1i = self.fc('1-f-i')
        f1o = self.fc('1-f-o')
        f2 = self.fc('2-f-o')

        w0 = self.fc('0-f-i')
        w1i = self.fc('1-f-i')
        w1o = self.fc('1-f-o')
        w2 = self.fc('2-f-o')

        M0 = f0.matrices.mass
        M2 = f2.matrices.mass
        M1i = f1i.matrices.mass
        M1o = f1o.matrices.mass

        W0 = f0.operators.inner(w0)
        W2 = f2.operators.inner(w2)
        W1i = f1i.operators.inner(w1i)
        W1o = f1o.operators.inner(w1o)

        for e in M0:
            np.testing.assert_array_almost_equal(M0[e].toarray(), W0[e].toarray())
            np.testing.assert_array_almost_equal(M2[e].toarray(), W2[e].toarray())
            np.testing.assert_array_almost_equal(M1i[e].toarray(), W1i[e].toarray())
            np.testing.assert_array_almost_equal(M1o[e].toarray(), W1o[e].toarray())

        return 1



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/unittests/standard_forms/mass_matrices.py
    miUsGrid_Triangles_MassMatrices()()
