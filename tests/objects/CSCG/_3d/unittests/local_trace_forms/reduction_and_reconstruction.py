# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 5:40 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

from components.miscellaneous.miprint import miprint

from __init__ import cscg3
import numpy as np


def p(t, x, y, z):
    return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.cos(2 * np.pi * z) + t


def u(t, x, y, z):
    return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.cos(2 * np.pi * z) + t


def v(t, x, y, z):
    return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.cos(2 * np.pi * z) + t


def w(t, x, y, z):
    return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.sin(2 * np.pi * z) + t


class Test_Reduction_and_Reconstruction_of_local_trace_forms(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        mesh = cscg3.mesh('crazy', c=0.0)([5, 6, 5])
        space = cscg3.space('polynomials')([5, 4, 5])
        self.fc = cscg3.form(mesh, space)
        miprint(">>> Test_Reduction_and_Reconstruction_of_local_trace_forms ...", flush=True)
        self._freeze_self_()

    def __call__(self):
        """"""
        scalar = self.fc('scalar', p)
        vector = self.fc('vector', [u, v, w])
        scalar.current_time = 0
        vector.current_time = 0

        ltf0 = self.fc('0-lt')
        ltf1 = self.fc('1-lt')
        ltf2 = self.fc('2-lt')

        ltf0.CF = scalar
        ltf1.CF = vector
        ltf2.CF = scalar

        ltf0.discretize()
        ltf2.discretize()

        xi = list((np.random.rand(7)-0.5) * 2)
        et = list((np.random.rand(7)-0.5) * 2)
        sg = list((np.random.rand(7)-0.5) * 2)
        xi.sort()
        et.sort()
        sg.sort()

        xyz0, v0 = ltf0.reconstruct(xi, et, sg)
        xyz2, v2 = ltf2.reconstruct(xi, et, sg)

        for i in xyz0:
            xyz = xyz0[i]
            val = v0[i]

            for side in xyz:
                x, y, z = xyz[side]
                value = val[side][0]
                assert np.max(np.abs(p(0, x, y, z) - value)) < 1e-4

            xyz = xyz2[i]
            val = v2[i]

            for side in xyz:
                x, y, z = xyz[side]
                value = val[side][0]
                assert np.max(np.abs(p(0, x, y, z) - value)) < 1e-3

        R0 = ltf0.do.make_reconstruction_matrix_on_grid(xi, et, sg)
        R2 = ltf2.do.make_reconstruction_matrix_on_grid(xi, et, sg)

        ll0 = ltf0.cochain.local_ESW
        ll2 = ltf2.cochain.local_ESW

        xyz0, v0 = ltf0.reconstruct(xi, et, sg, ravel=True)
        xyz2, v2 = ltf2.reconstruct(xi, et, sg, ravel=True)

        for i in R0:
            for side in 'NSWEBF':

                V = R0[i][side] @ ll0[i][side]

                np.testing.assert_array_almost_equal(V, v0[i][side][0])

                V = R2[i][side] @ ll2[i][side]

                np.testing.assert_array_almost_equal(V, v2[i][side][0])

        for ltf in (self.fc('0-lt', hybrid=False), self.fc('2-lt', hybrid=False)):

            ltf.CF = scalar

            ltf.discretize()

            xi = list((np.random.rand(7)-0.5) * 2)
            et = list((np.random.rand(7)-0.5) * 2)
            sg = list((np.random.rand(7)-0.5) * 2)

            xi.sort()
            et.sort()
            sg.sort()

            R0 = ltf.do.make_reconstruction_matrix_on_grid(xi, et, sg)

            ll0 = ltf.cochain.local_ESW

            v0 = ltf.reconstruct(xi, et, sg, ravel=True)[1]

            for i in R0:
                for side in 'NSWEBF':

                    V = R0[i][side] @ ll0[i][side]

                    np.testing.assert_array_almost_equal(V, v0[i][side][0])

        return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/local_trace_forms/reduction_and_reconstruction.py
    Test_Reduction_and_Reconstruction_of_local_trace_forms()()
