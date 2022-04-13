
import sys
if './' not in sys.path: sys.path.append('./')
import random


from objects.CSCG._3d.__init__ import mesh as mesh3
from objects.CSCG._3d.__init__ import space as space3
from objects.CSCG._3d.__init__ import form as form3
from objects.CSCG._3d.__init__ import exact_solution as es3

from root.config.main import *
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate
from tools.linear_algebra.linear_system.main import LinearSystem

def test_applying_strong_BC_for_Poisson_problem_NT():
    """"""
    mesh = mesh3('crazy',
                 c=0, bounds=([0.125, 1.125], [0.125, 1.125], [0.125, 1.125]))(
                        [5,4,3], EDM=None)
    space = space3('polynomials')([('Lobatto', 2), ('Lobatto', 3), ('Lobatto', 4)])
    FC = form3(mesh, space)
    ES = es3(mesh)('Poisson:sincos1')

    all_boundaries = mesh.boundaries.names
    if rAnk == mAster_rank:
        rn = random.randint(1,5)
        boundaries = random.sample(all_boundaries, rn)
    else:
        boundaries = None
    boundaries = cOmm.bcast(boundaries, root=mAster_rank)

    u_boundaries = boundaries
    if rAnk == mAster_rank:
        print(f"NTP [Applying_NT_BC_for_Poisson] hdMSEM@3d {u_boundaries}. ", flush=True)
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

    f.TW.func.body = ES.status.source_term
    f.TW.do.push_all_to_instant(0)
    f.discretize()
    B1 = f.cochain.EWC
    B1.gathering_matrix = f

    B2 = EWC_ColumnVector(mesh, t)
    B2.gathering_matrix = t

    b = concatenate([B0, B1, B2])

    LS = LinearSystem(A, b)

    u.TW.BC.body = ES.status.velocity
    u.TW.do.push_BC_to_instant(0)
    u.BC.valid_boundaries = u_boundaries
    t.BC.valid_boundaries = u_boundaries
    upc = u.BC.partial_cochain
    tpd = t.BC.partial_dofs

    LS.customize.apply_strong_BC(2, 0, tpd, upc)  # this is the test!!!!

    t.prime.TW.BC.body = ES.status.potential
    t.prime.TW.do.push_BC_to_instant(0)
    t.BC.valid_boundaries = p_boundaries
    tpc = t.BC.partial_cochain

    LS.customize.apply_strong_BC(2, 2, tpc)  # this is the test!!!!

    results = LS.solve('direct')()[0]
    results.do.distributed_to(u, p, t)

    u.TW.func.body = ES.status.velocity
    u.TW.do.push_all_to_instant(0)
    u_error_L2 = u.error.L()

    p.prime.TW.func.body = ES.status.potential
    p.prime.TW.do.push_all_to_instant(0)
    p_error_L2 = p.error.L()
    p_error_dH1 = p.error.dH(t, ES.status.velocity, 0)

    assert u_error_L2 < 0.2 and p_error_L2 < 0.021 and p_error_dH1 < 0.2

    return 1

if __name__ == '__main__':
    # mpiexec -n 4 python __tests__\unittests\linear_system\strong_BC_of_Poisson_non_trivial_BC.py

    test_applying_strong_BC_for_Poisson_problem_NT()
