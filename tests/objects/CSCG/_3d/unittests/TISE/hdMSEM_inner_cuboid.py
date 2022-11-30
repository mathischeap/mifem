# -*- coding: utf-8 -*-
"""
Here we use the hdMSEM to solve the inner-orientated version of the time independent Schrödinger
equation. We do this to test the hybridization of 0-forms.

"""

import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.__init__ import mesh as mesh3
from objects.CSCG._3d.__init__ import space as space3
from objects.CSCG._3d.__init__ import form as form3
from objects.CSCG._3d.__init__ import exact_solution as es3

from root.config.main import *
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.miLinearAlgebra.linearSystem.main import LinearSystem

from objects.CSCG.tools.distribute_local_cochain import distribute_local_cochain



def test_hdMSEM_Schrodinger_Inner_Cuboid():
    """"""
    mesh = mesh3('cuboid', region_layout=[3,2,2])([2, 3, 3], EDM=None)
    space = space3('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
    FC = form3(mesh, space)
    ES = es3(mesh)('TISE:sincos1', m=1e-66, E=1)

    all_boundaries = mesh.boundaries.names

    u_boundaries = ['North', 'East', 'Back']
    if RANK == MASTER_RANK:
        print(f"inSC [test_hdMSEM_Schrodinger_Inner_Cuboid] @ u_boundaries = {u_boundaries}. ", flush=True)
    p_boundaries = list()
    for b in all_boundaries:
        if b not in u_boundaries:
            p_boundaries.append(b)

    u = FC('1-f', is_hybrid = True)
    p = FC('0-f', is_hybrid = True)
    t = FC('0-adt')
    e = FC('0-e')

    V = ES.V(0, 0, 0, 0)

    E = ES.E
    alpha = ES._alpha_

    p.BC.CF = ES.wave_function
    p.BC.CF.current_time = 0
    p.BC.boundaries = p_boundaries

    t.prime.BC.CF = ES.flux.flux
    t.prime.BC.CF.current_time = 0
    t.BC.boundaries = u_boundaries

    I = u.matrices.identity
    E10 = p.matrices.incidence
    E01 = E10.T
    M0 = p.matrices.mass
    M1 = u.matrices.mass
    T0T = t.matrices.trace.T
    T, D, C, b2, eGM = p.special.hybrid_pairing(t, e)

    A = bmat([(        I,             -E10, None, None),
              (-E01 @ M1, ((E-V)/alpha)*M0,  T0T, None),
              (     None,                T,    D,    C),
              (     None,             None,  C.T, None)])

    A.gathering_matrices = [(u, p, t, eGM), (u, p, t, eGM)]

    b0 = EWC_ColumnVector(mesh, u)
    b0.gathering_matrix = u

    b1 = EWC_ColumnVector(mesh, p)
    b1.gathering_matrix = p

    b2.gathering_matrix = t

    b3 = EWC_ColumnVector(mesh, e)
    b3.gathering_matrix = eGM

    b = concatenate([b0, b1, b2, b3])

    LS = LinearSystem(A, b)

    results = LS.solve('Schur', rank=2, blocks=2)()[0]
    distribute_local_cochain(results, [u, p, t, e])

    p.CF = ES.wave_function
    p.CF.current_time = 0
    p_error_L2 = p.error.L()

    u.CF = ES.flux
    u.CF.current_time = 0
    u_error_L2 = u.error.L()

    du = FC('1-adf', u)
    u_error_dH1 = du.error.dH(t, ES.source_term)

    assert p_error_L2 < 0.006
    assert u_error_L2 < 0.17
    assert u_error_dH1 < 0.7

    return 1






if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/TISE/hdMSEM_inner_cuboid.py

    # for _ in range(100):
    test_hdMSEM_Schrodinger_Inner_Cuboid()