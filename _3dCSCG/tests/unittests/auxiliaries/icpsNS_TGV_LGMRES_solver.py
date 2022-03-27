# -*- coding: utf-8 -*-
"""
mpiexec -n 12 python _3dCSCG/TESTS/__unittest_scripts__/icpsNS_TGV_LGMRES_solver.py

This one actually is mainly used to test the LGMRES solver.

"""


import sys
if './' not in sys.path: sys.path.append('../__unittest_scripts__/')

from numpy import pi
from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.data_structures.global_matrix.main import LocallyFullVector
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import bmat, concatenate
from tools.linear_algebra.solvers.parallel.allocator import ParallelSolverDistributor
from tools.iterators.simple import SimpleIterator
from time import time
from root.config.main import *
from root.mifem.save import save
from screws.miscellaneous.timer import check_multiple_close, check_almost_in_range



# import warnings

def NoHy_TGV_NEW_LGMRES(N=2, k=4, t=15, steps=480, Re=500,
    atol=1e-5, m=50, K=3, maxiter=10, show_info=True, save_uw=True):
    """

    :param N:
    :param k:
    :param t:
    :param steps: How many time steps in total?
    :param Re:
    :param atol:
    :param m:
    :param K:
    :param maxiter:
    :param show_info:
    :param save_uw:
    :return:
    """

    t0 = 0
    dt = (t - t0) / steps

    L = 1
    V0 = 1
    nu = 0 if Re > 9999 else V0 * L / Re

    Re, N, k, steps = int(Re), int(N), int(k), int(steps)

    quad_degree = [2 * N, 2 * N, 2 * N]
    RDF_filename = 'icpNS_NH_RDF_LGMRES__TGV_NEW_Re{}_N{}K{}t{}Steps{}'.format(Re, N, k, t, steps)
    ITR_name     = 'icpNS_NH_ITR_LGMRES__TGV_NEW_Re{}_N{}K{}t{}Steps{}'.format(Re, N, k, t, steps)
    auto_save_frequency = 5
    monitor_factor = 1

    SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                        auto_save_frequency=auto_save_frequency,
                        monitor_factor=monitor_factor,
                        RDF_filename=RDF_filename,
                        name=ITR_name)

    mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=[(-pi,pi),(-pi,pi),(-pi,pi)])([k, k, k], show_info=show_info)
    space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N), ('Lobatto', N)], show_info=show_info)
    FC = FormCaller(mesh, space)
    Volume = mesh.domain.volume
    es = ExactSolutionSelector(mesh)('icpsNS:TGV1', nu=nu, L=L, V0=V0, show_info=show_info)
    P0 = FC('0-f', is_hybrid=False, orientation='inner', name='inner-total-pressure')
    u1 = FC('1-f', is_hybrid=False, orientation='inner', name='inner-velocity')
    w2 = FC('2-f', is_hybrid=False, orientation='inner', name='inner-vorticity')
    w1 = FC('1-f', is_hybrid=False, orientation='outer', name='outer-vorticity')
    u2 = FC('2-f', is_hybrid=False, orientation='outer', name='outer-velocity')
    P3 = FC('3-f', is_hybrid=False, orientation='outer', name='outer-total-pressure')
    u1.TW.func.do.set_func_body_as(es.status.velocity)
    u2.TW.func.do.set_func_body_as(es.status.velocity)
    w1.TW.func.do.set_func_body_as(es.status.vorticity)
    w2.TW.func.do.set_func_body_as(es.status.vorticity)

    # define used fundamental matrices ......
    M1 = u1.matrices.mass
    M2 = u2.matrices.mass
    M3 = P3.matrices.mass
    E10 = P0.coboundary.incidence_matrix
    E21 = u1.coboundary.incidence_matrix
    E32 = u2.coboundary.incidence_matrix
    E12 = E21.T
    E23 = E32.T
    CP1 = w1.special.cross_product_1f__ip_1f(u1, u1, quad_degree=quad_degree)
    CP2 = w2.special.cross_product_2f__ip_2f(u2, u2, quad_degree=quad_degree)

    # ... compute t0 co-chains and conditions ......
    u1.TW.current_time = t0
    u2.TW.current_time = t0
    w1.TW.current_time = t0
    w2.TW.current_time = t0
    u1.TW.do.push_all_to_instant()
    u2.TW.do.push_all_to_instant()
    w1.TW.do.push_all_to_instant()
    w2.TW.do.push_all_to_instant()
    u1.discretize()
    u2.discretize()
    w1.discretize()
    w2.discretize()
    u1.error.L()
    u2.error.L()
    w1.error.L()
    w2.error.L()
    KE1_t0 = 0.5 * u1.do.compute_L2_energy_with(M=M1) / Volume
    KE2_t0 = 0.5 * u2.do.compute_L2_energy_with(M=M2) / Volume
    H1_t0 = u1.do.compute_L2_energy_with(w1, M=M1)
    H2_t0 = u2.do.compute_L2_energy_with(w2, M=M2)
    E1_t0 = 0.5 * w1.do.compute_L2_energy_with(M=M1) / Volume
    E2_t0 = 0.5 * w2.do.compute_L2_energy_with(M=M2) / Volume
    du2 = u2.coboundary()
    du2.TW.func.do.set_func_body_as(es.status.divergence_of_velocity)
    du2.TW.current_time = t0
    du2.TW.do.push_all_to_instant()
    DIV_L2_error_t0 = du2.error.L()
    u1u2_diff_t0 = u2.do.compute_Ln_diff_from(u1)
    w1w2_diff_t0 = w2.do.compute_Ln_diff_from(w1)

    if save_uw:
        save([u1, u2, w1, w2], f'icpNS_NH_UUWW_TGV_NEW_Re{Re}_N{N}k{k}t{t}Steps{steps}_t0')

    # set up inner half integer time step systems ......
    E12M2E21 = E12 @ M2 @ E21
    lhs00_Hf = 2*M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21
    lhs00_Hf.gathering_matrices = (u1, u1)
    M1E10 = M1 @ E10
    M1E10.gathering_matrices = (u1, P0)
    E01M1 = M1E10.T
    iA = bmat([[ lhs00_Hf, M1E10 ],
                [-E01M1   , None ]])
    iA.customize.identify_global_row(-1)

    B0 = (2 * M1 / dt - 0.5 * CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC
    B0.gathering_matrix = u1
    B1 = EWC_ColumnVector(mesh, P0.num.basis)
    B1.gathering_matrix = P0
    ib = concatenate([B0, B1])

    X0_0 = u1.cochain.EWC
    X0_1 = EWC_ColumnVector(mesh, P0.num.basis)
    X0_1.gathering_matrix = P0
    X0 = concatenate((X0_0, X0_1))

    iA = iA.assembled
    ib = ib.assembled
    X0 = LocallyFullVector(X0.assembled)

    iR = ParallelSolverDistributor("LGMRES")(iA, ib, X0, atol=atol, m=m, k=K, maxiter=maxiter)[0]
    iR.___PRIVATE_be_distributed_to___(u1, P0)

    w2.cochain.local = u1.coboundary.cochain_local
    KE1_t0h = 0.5 * u1.do.compute_L2_energy_with(M=M1) / Volume
    E2_t0h  = 0.5 * w2.do.compute_L2_energy_with(M=M2) / Volume

    if rAnk == mAster_rank and show_info:
        print('KE1_t0', KE1_t0)
        print('KE2_t0', KE2_t0)
        print('E1_t0', E1_t0)
        print('E2_t0', E2_t0)
        print('H1_t0', H1_t0)
        print('H2_t0', H2_t0)
        print('KE1_t0h', KE1_t0h)
        print('E2_t0h', E2_t0h, flush=True)

    iA00 =  M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21
    iA00.gathering_matrices = (u1, u1)
    iB_0 = (M1/dt - 0.5*CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC
    iB_0.gathering_matrix = u1

    iA = bmat( [[ iA00   , M1E10 ],
                [-E01M1  , None ]])
    iA.customize.identify_global_row(-1)
    ib = concatenate([iB_0, B1])

    # set up outer integer time step systems ......
    oA00 = M2 / dt + 0.5 * CP2
    oA00.gathering_matrices = (u2, u2)
    E23M3 = E23 @ M3
    E23M3.gathering_matrices = (u2, P3)
    M2E21 = M2 @ E21
    M2E21.gathering_matrices = (u2, w1)
    E12M2 = E12 @ M2
    E12M2.gathering_matrices = (w1, u2)
    M1.gathering_matrices = (w1, w1)
    E32.gathering_matrices = (P3, u2)

    lhs = [[ oA00  , 0.5*nu*M2E21, -E23M3  ],  # u2
           [-E12M2 , M1          , None    ],  # w1
           [ E32   , None        , None    ]]  # P3


    # del E23M3, mE23M3_A, E12M2, mE12M2_A, M1_A, E32_A, M2E21_A, E12M2E21

    oA = bmat(lhs)

    B0 = (M2 / dt - 0.5 * CP2) @ u2.cochain.EWC - 0.5*nu*M2E21 @ w1.cochain.EWC
    B0.gathering_matrix = u2

    B1 = EWC_ColumnVector(mesh, w1.num.basis)
    B1.gathering_matrix = w1


    B2 = EWC_ColumnVector(mesh, P3.num.basis)
    B2.gathering_matrix = P3
    ob = concatenate([B0, B1, B2])

    OUT_R = [0, ]
    INN_R = [iR,]



    def SOLVER(tk, tk1):
        """
        Parameters
        ----------
        tk :
        tk1 :

        Returns
        -------
        exit_code: The standard exit code.
        shut_down: If it is ``True``, the outer iterator will shutdown immediately.
        message: The solver message.
        KE1_tk :
        KE1_tkh :
        KE2_tk :
        H1_tk :
        H2_tk :
        E1_tk :
        E2_tk :
        E2_tkh :
        u2u1_L2_diff :
        w2w1_L2_diff :
        DIV_L2_error :
        """
        ts = time()
        assert tk1 == tk + dt


        SYS_oA = oA.assembled
        SYS_ob = ob.assembled


        if tk == t0:  # first step
            X0_0 = u2.cochain.EWC
            X0_1 = w1.cochain.EWC
            X0_2 = EWC_ColumnVector(mesh, P3.num.basis)
            X0_2.gathering_matrix = P3
            X0 = LocallyFullVector(concatenate((X0_0, X0_1, X0_2)).assembled)
        else:
            X0 = OUT_R[0]

        oR,_,_,_,mo = ParallelSolverDistributor("LGMRES")(SYS_oA, SYS_ob, X0, atol=atol, m=m, k=K, maxiter=maxiter)

        del SYS_oA, SYS_ob

        OUT_R[0] = oR
        oR.___PRIVATE_be_distributed_to___(u2, w1, P3)

        du2 = u2.coboundary()
        du2.TW.func.do.set_func_body_as(es.status.divergence_of_velocity)
        du2.TW.current_time = tk1
        du2.TW.do.push_all_to_instant()
        DIV_L2_error_tk1 = du2.error.L()
        KE2_tk1 = 0.5 * u2.do.compute_L2_energy_with(M=M2) / Volume
        E1_tk1 = 0.5 * w1.do.compute_L2_energy_with(M=M1) / Volume

        # ... inner

        SYS_iA = iA.assembled
        SYS_ib = ib.assembled

        X0 = INN_R[0]
        iR,_,_,_,mi = ParallelSolverDistributor("LGMRES")(SYS_iA, SYS_ib, X0, atol=atol, m=m, k=K, maxiter=maxiter)

        del SYS_iA, SYS_ib

        INN_R[0] = iR
        _u1_old_cochain_ = u1.cochain.local
        iR.___PRIVATE_be_distributed_to___(u1, P0)

        _u1_new_cochain_ = u1.cochain.local

        mean_u1_cochain_local_at_tk = dict()
        for i in _u1_old_cochain_:
            mean_u1_cochain_local_at_tk[i] = (_u1_old_cochain_[i] + _u1_new_cochain_[i]) / 2

        u1.cochain.local = mean_u1_cochain_local_at_tk  # we then have u1 cochain @ tk
        KE1_tk1 = 0.5 * u1.do.compute_L2_energy_with(M=M1) / Volume
        H1_tk1 = u1.do.compute_L2_energy_with(w1, M=M1)
        u1u2_diff_tk1 = u2.do.compute_Ln_diff_from(u1)

        w2.cochain.local = u1.coboundary.cochain_local
        H2_tk1 = u2.do.compute_L2_energy_with(w2, M=M2)
        E2_tk1 = 0.5 * w2.do.compute_L2_energy_with(M=M2) / Volume
        w1w2_diff_tk1 = w2.do.compute_Ln_diff_from(w1)

        if save_uw:
            if check_multiple_close(tk1, 0.1) and check_almost_in_range(tk1, 8.8, 9.5):
                TK1 = round(tk1, 1)
                save([u1, u2, w1, w2], f'icpNS_NH_UUWW_TGV_NEW_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            elif check_multiple_close(tk1, 1):
                TK1 = round(tk1)
                save([u1, u2, w1, w2], f'icpNS_NH_UUWW_TGV_NEW_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            else:
                pass
        u1.cochain.local = _u1_new_cochain_  # renew u1 cochain to time tk+half
        w2.cochain.local = u1.coboundary.cochain_local  # renew w2 cochain to time tk+half
        KE1_tk1h = 0.5 * u1.do.compute_L2_energy_with(M=M1) / Volume
        E2_tk1h = 0.5 * w2.do.compute_L2_energy_with(M=M2) / Volume
        if save_uw:
            if check_multiple_close(tk1, 0.1) and check_almost_in_range(tk1, 8.8, 9.5):
                TK1 = round(tk1, 1)
                save([u1, w2], f'icpNS_NH_UWih_TGV_NEW_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            elif check_multiple_close(tk1, 1):
                TK1 = round(tk1)
                save([u1, w2], f'icpNS_NH_UWih_TGV_NEW_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            else:
                pass

        message = [f'ITERATION cost: {int((time() - ts) * 100) / 100}',
                   'Inner solver: ' + mi,
                   'Outer solver: ' + mo]

        # print(KE1_tk1, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1)
        # print(message)

        return 1, 0, message, KE1_tk1, KE1_tk1h, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, E2_tk1h, \
               u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1



    SI(SOLVER, [KE1_t0, KE1_t0h, KE2_t0, H1_t0, H2_t0, E1_t0, E2_t0, E2_t0h, u1u2_diff_t0,
                w1w2_diff_t0, DIV_L2_error_t0])
    SI.run()

    return SI

if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG/TESTS/__unittest_scripts__/icpsNS_TGV_LGMRES_solver.py

    NoHy_TGV_NEW_LGMRES(N=2, k=2, t=1, steps=5, Re=500, atol=1e-5, m=30, K=3, maxiter=3, save_uw=False)