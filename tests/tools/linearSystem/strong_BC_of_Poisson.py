# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
import random

from objects.CSCG._2d.__init__ import mesh as mesh2
from objects.CSCG._2d.__init__ import space as space2
from objects.CSCG._2d.__init__ import form as form2
from objects.CSCG._2d.__init__ import exact_solution as es2

from objects.CSCG._3d.__init__ import mesh as mesh3
from objects.CSCG._3d.__init__ import space as space3
from objects.CSCG._3d.__init__ import form as form3
from objects.CSCG._3d.__init__ import exact_solution as es3

from root.config.main import *
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from tools.elementwiseCache.dataStructures.operators.bmat.main import bmat
from tools.elementwiseCache.dataStructures.operators.concatenate.main import concatenate
from tools.miLinearAlgebra.linearSystem.main import LinearSystem


def test_applying_strong_BC_for_Poisson_problem():
    """"""
    #---------- 2d test --------------------------------------------
    mesh = mesh2('rectangle', region_layout=(2,2))([4, 5], EDM=None)
    space = space2('polynomials')([('Lobatto', 4), ('Lobatto', 3)])
    FC = form2(mesh, space)
    ES = es2(mesh)('sL:sincos1')

    all_boundaries = mesh.boundaries.names
    if RANK == MASTER_RANK:
        rn = random.randint(0,3)
        boundaries = random.sample(all_boundaries, rn)
        print(f"SBC [Applying_strong_BC_for_Poisson] MSEM@2d {boundaries}. ", flush=True)
    else:
        boundaries = None
    boundaries = COMM.bcast(boundaries, root=MASTER_RANK)

    u = FC('1-f-o', is_hybrid=False)
    p = FC('2-f-o', is_hybrid=False)
    f = FC('2-f-o', is_hybrid=False)

    B0 = EWC_ColumnVector(mesh, u)
    B0.gathering_matrix = u
    f.CF = ES.source_term
    f.CF.current_time = 0
    f.discretize()
    B1 = f.cochain.EWC
    B1.gathering_matrix = f

    M1 = u.matrices.mass
    M2 = p.matrices.mass
    E21 = u.matrices.incidence
    E12 = E21.T
    E12M2 = E12 @ M2

    A = bmat(([M1 , E12M2],
              [E21, None ]))
    A.gathering_matrices = ((u, p), (u, p))
    b = concatenate([B0, B1])
    LS = LinearSystem(A, b)

    u.BC.CF = ES.velocity
    u.BC.CF.current_time = 0
    u.BC.boundaries = boundaries

    LS.customize.apply_strong_BC(0, 0, u) # this is the test!!!!

    results = LS.solve('direct')()[0]
    results.do.distributed_to(u, p)

    u.CF = ES.velocity
    u.CF.current_time = 0
    u_error_L2 = u.error.L()
    p.CF = ES.potential
    p.CF.current_time = 0
    p_error_L2 = p.error.L()

    assert u_error_L2 < 0.003 and p_error_L2 < 0.0004

    #---------- 3d test --------------------------------------------
    mesh = mesh3('crazy', c=0)([5,4,3], EDM=None)
    space = space3('polynomials')([('Lobatto', 2), ('Lobatto', 3), ('Lobatto', 4)])
    FC = form3(mesh, space)
    ES = es3(mesh)('Poisson:sincos1')

    all_boundaries = mesh.boundaries.names
    if RANK == MASTER_RANK:
        rn = random.randint(1,5)
        boundaries = random.sample(all_boundaries, rn)
        print(f"    [Applying_strong_BC_for_Poisson] MSEM@3d {boundaries}. ", flush=True)
    else:
        boundaries = None
    boundaries = COMM.bcast(boundaries, root=MASTER_RANK)
    u = FC('2-f', is_hybrid=False)
    p = FC('3-f', is_hybrid=False)
    f = FC('3-f', is_hybrid=False)

    B0 = EWC_ColumnVector(mesh, u)
    B0.gathering_matrix = u
    f.CF = ES.source_term
    f.CF.current_time = 0
    f.discretize()
    B1 = f.cochain.EWC
    B1.gathering_matrix = f

    M2 = u.matrices.mass
    M3 = p.matrices.mass
    E32 = u.matrices.incidence
    E23 = E32.T
    E23M3 = E23 @ M3

    A = bmat(([  M2, E23M3],
              [-E32, None ]))
    A.gathering_matrices = ((u, p), (u, p))
    b = concatenate([B0, B1])
    LS = LinearSystem(A, b)

    u.BC.CF = ES.velocity
    u.BC.CF.current_time = 0
    u.BC.boundaries = boundaries

    LS.customize.apply_strong_BC(0, 0, u) # this is the test!!!!

    results = LS.solve('direct')()[0]
    results.do.distributed_to(u, p)

    u.CF = ES.velocity
    u.CF.current_time = 0
    u_error_L2 = u.error.L()
    p.CF = ES.potential
    p.CF.current_time = 0
    p_error_L2 = p.error.L()

    assert u_error_L2 < 0.2 and p_error_L2 < 0.022

    # ---------- 3d test hdMSEM --------------------------------------------
    u_boundaries = boundaries
    if RANK == MASTER_RANK:
        print(f"    [Applying_strong_BC_for_Poisson] hdMSEM@3d {u_boundaries}. ", flush=True)
    p_boundaries = list()
    for b in all_boundaries:
        if b not in u_boundaries:
            p_boundaries.append(b)

    u = FC('2-f', is_hybrid = True)
    p = FC('3-adf')
    f = FC('3-f', is_hybrid = True)
    t = FC('2-adt')

    M2 = u.matrices.mass
    E32 = u.matrices.incidence
    E23 = E32.T
    T = t.coboundary.trace_matrix

    A = bmat(([  M2, E23 , -T.T],
              [-E32, None, None],
              [ T  , None, None]))
    A.gathering_matrices = ((u, p, t), (u, p, t))

    B0 = EWC_ColumnVector(mesh, u)
    B0.gathering_matrix = u

    f.CF = ES.source_term
    f.CF.current_time = 0
    f.discretize()
    B1 = f.cochain.EWC
    B1.gathering_matrix = f

    B2 = EWC_ColumnVector(mesh, t)
    B2.gathering_matrix = t

    b = concatenate([B0, B1, B2])

    LS = LinearSystem(A, b)

    u.BC.CF = ES.velocity
    u.BC.CF.current_time = 0
    u.BC.boundaries = u_boundaries
    t.BC.boundaries = u_boundaries
    upc = u.BC.interpret
    tpd = t.BC.interpret

    LS.customize.apply_strong_BC(2, 0, tpd, upc)  # this is the test!!!!

    t.prime.BC.CF = ES.potential
    ES.potential.current_time = 0
    t.BC.boundaries = p_boundaries
    tpc = t.BC.interpret

    LS.customize.apply_strong_BC(2, 2, tpc)  # this is the test!!!!

    results = LS.solve('direct')()[0]
    results.do.distributed_to(u, p, t)

    u.CF = ES.velocity
    u.CF.current_time = 0
    u_error_L2 = u.error.L()

    p.prime.CF = ES.potential
    p.prime.CF.current_time = 0
    p_error_L2 = p.error.L()
    p_error_dH1 = p.error.dH(t, ES.velocity, 0)

    assert u_error_L2 < 0.2 and p_error_L2 < 0.022 and p_error_dH1 < 0.2

    return 1






if __name__ == '__main__':
    # mpiexec -n 4 python __tests__/unittests/linear_system/strong_BC_of_Poisson.py

    test_applying_strong_BC_for_Poisson_problem()
