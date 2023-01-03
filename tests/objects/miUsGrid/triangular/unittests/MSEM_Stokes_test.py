# -*- coding: utf-8 -*-
"""
w = curl u

curl w + grad p = f

div u = 0

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/10/2022 1:26 AM
"""

import sys
if './' not in sys.path:
    sys.path.append('./')

from __init__ import miTri
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.miLinearAlgebra.linearSystem.main import LinearSystem
from root.config.main import RANK, MASTER_RANK

import numpy as np


def _Stokes2D_Solver_miUsTriangle(K=16, N=3):

    FC = miTri.call(f'st{K}', N)

    es = FC('Stokes : sin cos 1')

    w = FC('0-f-o', name='vorticity')
    u = FC('1-f-o', name='velocity')
    p = FC('2-f-o', name='pressure')
    f = FC('1-f-o', name='body_force')

    w.CF = es.vorticity
    u.CF = es.velocity
    p.CF = es.pressure
    f.CF = es.body_force

    u.BC.CF = es.velocity

    u.BC.boundaries = []

    es.current_time = 0  # set current time to all fields.

    f.discretize()

    M1 = w.matrices.mass
    M2 = u.matrices.mass
    M3 = p.matrices.mass
    E21 = w.matrices.incidence
    E32 = u.matrices.incidence
    E12 = E21.T
    E23 = E32.T

    A = ([M1,       - E12 @ M2,  None],
         [M2 @ E21, None,       -E23 @ M3],
         [None,     M3 @ E32,    None],)
    A = bmat(A)

    b0 = EWC_ColumnVector(w.mesh.elements, w.num.basis)
    b1 = M2 @ f
    b2 = EWC_ColumnVector(p.mesh.elements, p.num.basis)
    b = concatenate([b0, b1, b2])

    A.assembler.chain_method = 'sequent'
    b.assembler.chain_method = 'sequent'

    A.gathering_matrices = ([w, u, p], [w, u, p])
    b.gathering_matrix = (w, u, p)

    LS = LinearSystem(A, b)
    LS.customize.apply_strong_BC(1, 1, u)

    results = LS.solve('direct')()[0]

    results.do.distributed_to(w, u, p, chain_method='sequent')

    w_L2 = w.error.L()
    u_L2 = u.error.L()
    p_L2 = p.error.L()

    return w_L2, u_L2, p_L2


def miUsTriangleTest_MSEM_STOKES():
    if RANK == MASTER_RANK:
        print(f"--- [miUsTriangleTest_MSEM_STOKES] ...", flush=True)

    w0, u0, p0 = _Stokes2D_Solver_miUsTriangle(K=8, N=3)
    w1, u1, p1 = _Stokes2D_Solver_miUsTriangle(K=10, N=3)

    denominator = np.log10(1/8) - np.log10(1/10)

    numerator_w = np.log10(w0) - np.log10(w1)
    numerator_u = np.log10(u0) - np.log10(u1)
    numerator_p = np.log10(p0) - np.log10(p1)

    order_w = numerator_w / denominator
    order_u = numerator_u / denominator
    order_p = numerator_p / denominator

    np.testing.assert_almost_equal(order_w, 4, decimal=1)
    np.testing.assert_almost_equal(order_u, 3, decimal=1)
    np.testing.assert_almost_equal(order_p, 2, decimal=1)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/unittests/MSEM_Stokes_test.py
    miUsTriangleTest_MSEM_STOKES()
