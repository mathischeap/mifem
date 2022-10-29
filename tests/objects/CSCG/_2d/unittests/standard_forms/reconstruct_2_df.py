# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/30 3:57 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np

from __init__ import screws
from __init__ import cscg2


def test_reconstruct_S0F_to_DV():
    """"""
    screws.miprint("S0f [test_reconstruct_S0F_to_DV] ...... ", flush=True)
    mesh = cscg2.mesh('crazy', bounds=[(-1,1) for _ in range(2)])([8,9], EDM='debug')
    space = cscg2.space('polynomials')([6,5])
    FC = cscg2.form(mesh, space)
    def U(t, x, y): return 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y)
    f0 = FC('0-f-o')
    f2 = FC('2-f-i')
    Sca = FC('scalar', U)

    f0.CF = Sca
    f0.CF.current_time = 1
    f0.do.discretize()
    L2norm_0 = f0.do.compute_Ln_energy(n=2)

    f2.CF = Sca
    f2.CF.current_time = 1
    f2.do.discretize()
    L2norm_2 = f2.do.compute_Ln_energy(n=2)

    Qn, Qw = screws.Quadrature([14,15]).quad

    DS = f0.reconstruct.discrete_scalar(Qn)
    coordinates = DS.coordinates
    values = DS.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ij, i, j ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=6)
    DS = f2.reconstruct.discrete_scalar(Qn)
    coordinates = DS.coordinates
    values = DS.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ij, i, j ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_2, L2_norm, decimal=5)

    Qn, Qw = screws.Quadrature([5,6]).quad

    DS = f0.reconstruct.discrete_scalar(Qn)
    coordinates = DS.coordinates
    values = DS.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ij, i, j ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_0, L2_norm, decimal=2)
    DS = f2.reconstruct.discrete_scalar(Qn)
    coordinates = DS.coordinates
    values = DS.values
    for rn in coordinates:
        value = values[rn]
        v0 = value[0]
        L2_norm = np.einsum('ij, i, j ->', v0**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_2, L2_norm, decimal=2)

    return 1

def test_reconstruct_oS1F_to_DV():
    """"""
    screws.miprint("o1f [test_reconstruct_oS1F_to_DV] ...... ", flush=True)
    mesh = cscg2.mesh('crazy', bounds=[(-1,1) for _ in range(2)])([7,8], EDM='debug')
    space = cscg2.space('polynomials')([5,4])
    FC = cscg2.form(mesh, space)
    def U(t, x, y): return 2 * t * np.sin(np.pi * x) * np.cos(np.pi * y)
    def V(t, x, y): return t * np.sin(2 * np.pi * x) * np.cos(np.pi * y)
    u = FC('1-f-o')
    v = FC('1-f-i')
    vector = FC('vector', [U, V])
    u.CF = vector
    u.CF.current_time = 1
    u.do.discretize()
    L2norm_u = u.do.compute_Ln_energy(n=2)
    v.CF = vector
    v.CF.current_time = 1
    v.do.discretize()
    L2norm_v = v.do.compute_Ln_energy(n=2)

    Qn, Qw = screws.Quadrature([14,15]).quad
    DV = u.reconstruct.discrete_vector(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0, v1 = value
        L2_norm = np.einsum('ij, i, j ->', v0**2 + v1**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_u, L2_norm, decimal=3)
    DV = v.reconstruct.discrete_vector(Qn)
    coordinates = DV.coordinates
    values = DV.values
    for rn in coordinates:
        value = values[rn]
        v0, v1 = value
        L2_norm = np.einsum('ij, i, j ->', v0**2 + v1**2, *Qw, optimize='greedy')
        np.testing.assert_almost_equal(L2norm_v, L2_norm, decimal=3)

    Qn, Qw = screws.Quadrature([6, 7]).quad
    DV_coarse = u.reconstruct.discrete_vector(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0, v1 = value
        L2_norm = np.einsum('ij, i, j ->', v0**2 + v1**2, *Qw, optimize='greedy')
        assert L2norm_u - L2_norm < 0.2
    DV_coarse = v.reconstruct.discrete_vector(Qn)
    coordinates = DV_coarse.coordinates
    values = DV_coarse.values
    for rn in coordinates:
        value = values[rn]
        v0, v1 = value
        L2_norm = np.einsum('ij, i, j ->', v0**2 + v1**2, *Qw, optimize='greedy')
        assert L2norm_v - L2_norm < 0.2

    return 1


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/standard_forms/reconstruct_2_df.py
    test_reconstruct_oS1F_to_DV()
