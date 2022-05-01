
import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat
from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector
from tools.linear_algebra.linear_system.main import LinearSystem
from root.config.main import rAnk, mAster_rank
from objects.CSCG.tools.distribute_local_cochain import distribute_local_cochain



def test_Stokes_hdMSEM_Schur_Rank2Solver():
    if rAnk == mAster_rank:
        print(f"STK [test_Stokes_hdMSEM_Schur_Rank2Solver] ...", flush=True)

    K = [4, 3, 2]
    N = [2, 3, 4]
    mesh = MeshGenerator('crazy', c=0.0,
        bounds=[(0.125, 1.125),(0.125, 1.125),(0.125, 1.125)])(
        K, show_info=True)

    space = SpaceInvoker('polynomials')(N, show_info=True)
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('Stokes:sincos1')

    p_Boundaries = ['Back', 'North', 'South', ]
    n_u_Boundaries = ['Front', 'West', "East", ]

    t_u_Boundaries = ['Back', 'North', 'South', ]
    w_Boundaries = ['Front', 'West', "East", ]

    w = FC('1-f', is_hybrid=True, name='vorticity')
    u = FC('2-f', is_hybrid=True, name='velocity')
    p = FC('3-adf', name='pressure')
    f = FC('2-f', is_hybrid=True, name='body_force')

    t = FC('2-adt', name='pressure_trace')
    s = FC('1-adt', name='tangential_velocity')
    e = FC('1-e', name='edge-form')

    w.TW.func.do.set_func_body_as(es, 'vorticity')
    u.TW.func.do.set_func_body_as(es, 'velocity')
    p.prime.TW.func.do.set_func_body_as(es, 'pressure')
    f.TW.func.do.set_func_body_as(es, 'body_force')

    T2T = t.matrices.trace.T
    T1T = s.matrices.trace.T

    w.TW.BC.body = es.status.vorticity
    w.BC.valid_boundaries = w_Boundaries
    s.prime.TW.BC.body = es.status.velocity.components.T_perp
    s.BC.valid_boundaries = t_u_Boundaries
    T1, D1, C1, b3, eGM = w.special.hybrid_pairing(s, e, time=0)

    u.TW.BC.body = es.status.velocity
    u.BC.valid_boundaries = n_u_Boundaries
    t.prime.TW.BC.body = es.status.pressure
    t.BC.valid_boundaries = p_Boundaries
    T2, D2, b4 = u.special.hybrid_pairing(t, time=0)

    M1 = w.matrices.mass
    M2 = u.matrices.mass
    E21 = w.matrices.incidence
    E32 = u.matrices.incidence
    E12 = E21.T
    E23 = E32.T

    A = ([      M1, - E12 @ M2,  None,  T1T, None, None],
         [M2 @ E21,       None, - E23, None,  T2T, None],
         [    None,        E32,  None, None, None, None],
         [      T1,       None,  None,   D1, None,   C1],
         [    None,         T2,  None, None,   D2, None],
         [    None,       None,  None, C1.T, None, None])

    A = bmat(A)

    f.TW.current_time = 0
    f.TW.do.push_func_to_instant()
    f.discretize()

    b0 = EWC_ColumnVector(w)
    b1 = M2 @ f
    b2 = EWC_ColumnVector(p)
    b5 = EWC_ColumnVector(e)

    b = concatenate([b0, b1, b2, b3, b4, b5])

    A.gathering_matrices = ([w, u, p, s, t, eGM], [w, u, p, s, t, eGM])
    b.gathering_matrix = [w, u, p, s, t, eGM]

    LS = LinearSystem(A, b)

    results = LS.solve('Schur', rank=2, blocks=3)()[0]

    distribute_local_cochain(results, [w, u, p, s, t, e])

    w.TW.current_time = 0
    w.TW.do.push_func_to_instant()
    u.TW.current_time = 0
    u.TW.do.push_func_to_instant()
    p.prime.TW.current_time = 0
    p.prime.TW.do.push_func_to_instant()
    f.TW.current_time = 0
    w_L2 = w.error.L()
    u_L2 = u.error.L()
    p_L2 = p.error.L()

    p_error_dH1 = p.error.dH(t, es.status.gradient_of_pressure)
    u_error_dH1 = u.error.dH(s, es.status.vorticity)
    u_error_H1 = u.error.H(es.status.divergence_of_velocity)

    du_L2 = u.do.compute_Ln_norm_of_coboundary()

    assert du_L2 < 1e-10, f"conservation of mass."
    assert w_L2 < 0.7
    assert u_L2 < 0.08
    assert u_error_H1 < 0.08
    assert p_L2 < 0.04
    assert p_error_dH1 < 0.4
    assert u_error_dH1 < 0.7

    return 1





if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/Stokes_flow/hdMSEM_Schur_direct_3.py
    test_Stokes_hdMSEM_Schur_Rank2Solver()