# -*- coding: utf-8 -*-
"""
Here we use the hdMSEM to solve the inner-orientated version of the Poisson problem. We do this to
test the hybridization of 0-forms.

"""

import sys
if './' not in sys.path: sys.path.append('./')
import random

from objects.CSCG._3d.__init__ import mesh as mesh3
from objects.CSCG._3d.__init__ import space as space3
from objects.CSCG._3d.__init__ import form as form3
from objects.CSCG._3d.__init__ import exact_solution as es3

from root.config.main import *
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.miLinearAlgebra.linearSystem.main import LinearSystem



def test_hdMSEM_Poisson_Inner():
    """"""
    mesh = mesh3('crazy', c=0, bounds=([0.125, 1.125], [0.125, 1.125], [0.125, 1.125]))(
        [4, 3, 5], EDM=None)
    space = space3('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = form3(mesh, space)
    ES = es3(mesh)('Poisson:sincos1')

    all_boundaries = mesh.boundaries.names
    if RANK == MASTER_RANK:
        rn = random.randint(1, 5)
        boundaries = random.sample(all_boundaries, rn)
    else:
        boundaries = None
    boundaries = COMM.bcast(boundaries, root=MASTER_RANK)

    u_boundaries = boundaries
    if RANK == MASTER_RANK:
        print(f"inP [hdMSEM_inner_Poisson] @ u_boundaries = {u_boundaries}. ", flush=True)
    p_boundaries = list()
    for b in all_boundaries:
        if b not in u_boundaries:
            p_boundaries.append(b)

    p = FC('0-f', hybrid = True)
    u = FC('1-f', hybrid = True)
    t = FC('0-adt')
    e = FC('0-e')
    f = FC('0-adf')

    p.BC.CF = ES.potential
    p.BC.CF.current_time = 0
    p.BC.boundaries = p_boundaries

    t.prime.BC.CF = ES.velocity.flux
    t.prime.BC.CF.current_time = 0
    t.BC.boundaries = u_boundaries

    I = EWC_SparseMatrix(mesh, ('identity', u.num.basis))
    E10 = p.matrices.incidence
    E01 = E10.T
    M1 = u.matrices.mass
    T0T = t.matrices.trace.T
    T, D, C, b2, eGM = p.special.hybrid_pairing(t, e)

    A = bmat([(        I, -E10, None, None),
              (-E01 @ M1, None,  T0T, None),
              (     None,    T,    D,    C),
              (     None, None,  C.T, None)])
    A.gathering_matrices = [(u, p, t, eGM), (u, p, t, eGM)]

    b0 = EWC_ColumnVector(mesh, u)
    b0.gathering_matrix = u

    f.prime.CF = ES.source_term
    f.prime.CF.current_time = 0
    f.prime.discretize()

    b1 = - f.cochain.EWC
    b1.gathering_matrix = p

    b2.gathering_matrix = t

    b3 = EWC_ColumnVector(mesh, e)
    b3.gathering_matrix = eGM

    b = concatenate([b0, b1, b2, b3])

    LS = LinearSystem(A, b)

    results = LS.solve('direct')()[0]
    results.do.distributed_to(u, p, t)

    p.CF = ES.potential
    p.CF.current_time = 0
    p_error_L2 = p.error.L()

    u.CF = ES.velocity
    u.CF.current_time = 0
    u_error_L2 = u.error.L()

    du = FC('1-adf', u)
    u_error_dH1 = du.error.dH(t, - ES.source_term)

    f_error_L2 = f.prime.error.L()

    assert p_error_L2 < 0.005
    assert u_error_L2 < 0.14
    assert u_error_dH1 < 0.5
    assert f_error_L2 < 0.5

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/Poisson/hdMSEM_inner.py

    test_hdMSEM_Poisson_Inner()