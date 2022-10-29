# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 12:27 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np

from __init__ import screws
from __init__ import cscg3


def test_reconstruct_DF():
    """"""
    screws.miprint("3DF [test_reconstruct_DF] ...... ", flush=True)
    mesh = cscg3.mesh('crazy', bounds=[(-1,1) for _ in range(3)])([7,8,6], EDM='debug')
    space = cscg3.space('polynomials')([('Lobatto', 5), ('Lobatto', 4), ('Lobatto', 6)])
    FC = cscg3.form(mesh, space)
    def U(t, x, y, z): return 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    def V(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
    def W(t, x, y, z): return t * np.sin(2 * np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)

    vector = FC('vector', [U, V, W])
    scalar = FC('scalar', U)

    #---------------------------- standard 2 form ------------------------------------------------
    u = FC('2-f')

    Qn, Qw = screws.Quadrature([14, 15, 16]).quad
    u.CF = vector
    u.CF.current_time = 1
    u.do.discretize()
    L2norm_0 = u.do.compute_Ln_energy(n=2)
    DV = u.reconstruct.discrete_vector(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0, v1, v2 = value
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2 + v1**2 + v2**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=3)

    Qn, Qw = screws.Quadrature([7, 5, 6]).quad
    DV_coarse = u.reconstruct.discrete_vector(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0, v1, v2 = value
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2 + v1**2 + v2**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=1)

    #---------------------------- standard 1 form ------------------------------------------------
    u = FC('1-f')

    Qn, Qw = screws.Quadrature([14, 15, 16]).quad
    u.CF = vector
    u.CF.current_time = 1
    u.do.discretize()
    L2norm_0 = u.do.compute_Ln_energy(n=2)
    DV = u.reconstruct.discrete_vector(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0, v1, v2 = value
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2 + v1**2 + v2**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=4)

    Qn, Qw = screws.Quadrature([7, 5, 6]).quad
    DV_coarse = u.reconstruct.discrete_vector(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0, v1, v2 = value
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2 + v1**2 + v2**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=1)

    #------------- Standard 0-form --------------------------------------------------------------
    u = FC('0-f')

    Qn, Qw = screws.Quadrature([14, 15, 16]).quad
    u.CF = scalar
    u.CF.current_time = 1
    u.do.discretize()
    L2norm_0 = u.do.compute_Ln_energy(n=2)
    DV = u.reconstruct.discrete_field(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=4)

    Qn, Qw = screws.Quadrature([7, 5, 6]).quad
    DV_coarse = u.reconstruct.discrete_field(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=2)

    #------------- Standard 2-form --------------------------------------------------------------
    u = FC('3-f')

    Qn, Qw = screws.Quadrature([14, 15, 16]).quad
    u.CF = scalar
    u.CF.current_time = 1
    u.do.discretize()
    L2norm_0 = u.do.compute_Ln_energy(n=2)
    DV = u.reconstruct.discrete_field(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=3)

    Qn, Qw = screws.Quadrature([7, 5, 6]).quad
    DV_coarse = u.reconstruct.discrete_field(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ijk, i, j, k ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=2)

    return 1


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/standard_forms/reconstruct_2_DF.py
    test_reconstruct_DF()
