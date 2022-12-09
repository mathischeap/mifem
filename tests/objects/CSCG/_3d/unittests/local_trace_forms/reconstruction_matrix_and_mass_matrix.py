# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/29/2022 5:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from components.miscellaneous.miprint import miprint

from __init__ import cscg3
import numpy as np

class Test_reconstruction_matrix_and_mass_matrix(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        mesh = cscg3.mesh('crazy', c=0.0)([3, 4, 2])
        space = cscg3.space('polynomials')([2, 1, 3])
        self.fc = cscg3.form(mesh, space)
        miprint(">>> Test_reconstruction_matrix_and_mass_matrix ...", flush=True)
        self._freeze_self_()

    def __call__(self):
        """"""
        ltf0 = self.fc('0-lt')
        ltf2 = self.fc('2-lt')

        M0 = ltf0.matrices.mass
        M2 = ltf2.matrices.mass

        m0 = ltf0.___PrLT_mass_matrices_brutal_force___()
        m2 = ltf2.___PrLT_mass_matrices_brutal_force___()

        for i in m0:
            m = m0[i]
            M = M0[i]
            np.testing.assert_array_almost_equal(m.toarray(), M.toarray())

            m = m2[i]
            M = M2[i]
            np.testing.assert_array_almost_equal(m.toarray(), M.toarray())

        ltf0 = self.fc('0-lt', hybrid=False)
        ltf2 = self.fc('2-lt', hybrid=False)

        M0 = ltf0.matrices.mass
        M2 = ltf2.matrices.mass

        m0 = ltf0.___PrLT_mass_matrices_brutal_force___()
        m2 = ltf2.___PrLT_mass_matrices_brutal_force___()

        for i in m0:
            m = m0[i]
            M = M0[i]
            np.testing.assert_array_almost_equal(m.toarray(), M.toarray())

            m = m2[i]
            M = M2[i]
            np.testing.assert_array_almost_equal(m.toarray(), M.toarray())

        return 1

if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/local_trace_forms/reconstruction_matrix_and_mass_matrix.py
    Test_reconstruction_matrix_and_mass_matrix()()
