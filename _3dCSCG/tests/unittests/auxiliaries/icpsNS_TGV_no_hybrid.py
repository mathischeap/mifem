# -*- coding: utf-8 -*-
"""
mpiexec python _3dCSCG\APP\contents\icpsNS\no_hybrid\TGV.py

"""


import sys
if './' not in sys.path: sys.path.append('../__unittest_scripts__/')

from numpy import pi
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.data_structures import GlobalMatrix, GlobalVector, DistributedVector
from scipy import sparse as spspa
from tools.iterator import SimpleIterator
import tools.linear_algebra.deprecated.operators as TLO
import tools.linear_algebra.solvers.serial.scipy_sparse_linalg as scipy_sparse_linalg
from time import time
from root.config import *
from root.mifem import save
from screws.miscellaneous import check_multiple_close, check_almost_in_range

# import warnings

def NoHy_TGV(N=2, k=4, t=15, steps=480, Re=500,
    tol=1e-5, restart=50, maxiter=10, solver='gcrotmk', show_info=True, save_uw=True):
    """

    :param N:
    :param k:
    :param t:
    :param steps:
    :param Re:
    :param tol:
    :param restart:
    :param maxiter:
    :param solver:
    :param show_info:
    :param save_uw:
    :return:
    """

    t0 = 0
    dt = (t - t0) / steps

    L = 1
    V0 = 1
    nu = 0 if Re > 9999 else V0 * L / Re
    rho = 1

    Re, N, k, steps = int(Re), int(N), int(k), int(steps)

    quad_degree = [2 * N, 2 * N, 2 * N]
    RDF_filename = 'RDF_TGV_Re{}_N{}K{}t{}Steps{}'.format(Re, N, k, t, steps)
    ITR_name     = 'ITR_TGV_Re{}_N{}K{}t{}Steps{}'.format(Re, N, k, t, steps)
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
    es = ExactSolutionSelector(mesh)('icpsNS:TGV1', nu=nu, rho=rho, L=L, V0=V0, show_info=show_info)
    P0 = FC('0-f', is_hybrid=False, orientation='inner', name='inner-total-pressure')
    u1 = FC('1-f', is_hybrid=False, orientation='inner', name='inner-velocity')
    w2 = FC('2-f', is_hybrid=False, orientation='inner', name='inner-vorticity')
    w1 = FC('1-f', is_hybrid=False, orientation='outer', name='outer-vorticity')
    u2 = FC('2-f', is_hybrid=False, orientation='outer', name='outer-velocity')
    P3 = FC('3-f', is_hybrid=False, orientation='outer', name='outer-total-pressure')
    u1.TW.func.___DO_set_func_body_as___(es.status.velocity)
    u2.TW.func.___DO_set_func_body_as___(es.status.velocity)
    w1.TW.func.___DO_set_func_body_as___(es.status.vorticity)
    w2.TW.func.___DO_set_func_body_as___(es.status.vorticity)

    # define used fundamental matrices ......
    M1 = u1.matrices.mass
    M2 = u2.matrices.mass
    M3 = P3.matrices.mass
    E10 = P0.coboundary.incidence_matrix
    E21 = u1.coboundary.incidence_matrix
    E32 = u2.coboundary.incidence_matrix
    E12 = E21.T
    E23 = E32.T
    CP1 = w1.special.cross_product(u1, u1, quad_degree=quad_degree)
    CP2 = w2.special.cross_product(u2, u2, quad_degree=quad_degree)

    # ... compute t0 co-chains and conditions ......
    u1.TW.current_time = t0
    u2.TW.current_time = t0
    w1.TW.current_time = t0
    w2.TW.current_time = t0
    u1.TW.___DO_push_all_to_instant___()
    u2.TW.___DO_push_all_to_instant___()
    w1.TW.___DO_push_all_to_instant___()
    w2.TW.___DO_push_all_to_instant___()
    u1.discretize()
    u2.discretize()
    w1.discretize()
    w2.discretize()
    u1.error.L()
    u2.error.L()
    w1.error.L()
    w2.error.L()
    KE1_t0 = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume
    KE2_t0 = 0.5 * u2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume
    H1_t0 = u1.DO.compute_L2_inner_product_energy_with(w1, M=M1)
    H2_t0 = u2.DO.compute_L2_inner_product_energy_with(w2, M=M2)
    E1_t0 = 0.5 * w1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume
    E2_t0 = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume
    du2 = u2.coboundary()
    du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
    du2.TW.current_time = t0
    du2.TW.___DO_push_all_to_instant___()
    DIV_L2_error_t0 = du2.error.L()
    u1u2_diff_t0 = u2.DO.compute_L2_diff_from(u1)
    w1w2_diff_t0 = w2.DO.compute_L2_diff_from(w1)

    if save_uw:
        save([u1, u2, w1, w2], f'UUWW_TGV_Re{Re}_N{N}k{k}t{t}Steps{steps}_t0')

    # set up inner half integer time step systems ......
    E12M2E21 = E12 @ M2 @ E21
    lhs00_Hf = 2*M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21
    lhs00_Hf.gathering_matrices = (u1, u1)
    lhs00_Hf_A = lhs00_Hf.assembled
    M1E10 = M1 @ E10
    M1E10.gathering_matrices = (u1, P0)
    M1E10_A = M1E10.assembled
    E01M1_A = M1E10_A.T
    E01M1_A._M_ = E01M1_A.M.tolil()
    E01M1_A.M[-1,:] = 0
    lhs11_A = GlobalMatrix((P0.GLOBAL_num_dofs, P0.GLOBAL_num_dofs))
    lhs11_A.M[-1,-1] = 1
    lhs = [[ lhs00_Hf_A, M1E10_A],
           [-E01M1_A   , lhs11_A]]
    del lhs00_Hf, lhs00_Hf_A, M1E10, M1E10_A, E01M1_A, lhs11_A
    iA = TLO.bmat(lhs, format='csr')
    del lhs
    iA = iA.___PRIVATE_gather_M_to_core___(clean_local=True)
    iA = GlobalMatrix(iA)
    assert iA.IS_master_dominating

    B0 = (2 * M1 / dt - 0.5 * CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC
    B0.gathering_matrix = u1
    B0 = B0.assembled
    B1 = GlobalVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
    ib = TLO.concatenate([B0, B1])
    ib = ib.___PRIVATE_gather_V_to_core___(clean_local=True)
    ib = GlobalVector(ib)
    assert ib.IS_master_dominating
    del B0, B1

    X0_0 = u1.cochain.globe
    X0_1 = DistributedVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
    X0 = TLO.concatenate((X0_0, X0_1))

    iR = getattr(scipy_sparse_linalg, solver)(
        iA, ib, X0, tol=tol, restart=restart, maxiter=maxiter)[0]
    iR.DO_distribute_to(u1, P0)
    w2.cochain.local = u1.coboundary.cochain_local
    KE1_t0h = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume
    E2_t0h  = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume

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

    # set up outer integer time step systems ......
    oA00 = M2 / dt + 0.5 * CP2
    oA00.gathering_matrices = (u2, u2)
    E23M3 = E23 @ M3
    E23M3.gathering_matrices = (u2, P3)
    mE23M3_A = - E23M3.assembled
    M2E21 = M2 @ E21
    M2E21.gathering_matrices = (u2, w1)
    M2E21_A = M2E21.assembled
    E12M2 = E12 @ M2
    E12M2.gathering_matrices = (w1, u2)
    mE12M2_A = - E12M2.assembled
    M1.gathering_matrices = (w1, w1)
    M1_A = M1.assembled
    E32.gathering_matrices = (P3, u2)
    E32_A = E32.assembled
    lhs = [[None    , 0.5*nu*M2E21_A, mE23M3_A  ],  # u2
           [mE12M2_A, M1_A          , None      ],  # w1
           [E32_A   , None          , None      ]]  # P3
    del E23M3, mE23M3_A, E12M2, mE12M2_A, M1_A, E32_A, M2E21_A, E12M2E21
    oA = TLO.bmat(lhs, format='csr')
    del lhs
    oA = oA.___PRIVATE_gather_M_to_core___(clean_local=True)
    oA = GlobalMatrix(oA)
    assert oA.IS_master_dominating

    oB_0 = (M2 / dt - 0.5 * CP2) @ u2.cochain.EWC - 0.5*nu*M2E21 @ w1.cochain.EWC
    oB_0.gathering_matrix = u2
    B0 = oB_0.assembled
    B1 = GlobalVector(spspa.csc_matrix((w1.GLOBAL_num_dofs, 1)))
    B2 = GlobalVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
    ob = TLO.concatenate([B0, B1, B2])
    ob = ob.___PRIVATE_gather_V_to_core___(clean_local=True)
    ob = GlobalVector(ob)
    assert ob.IS_master_dominating
    del B0, B1, B2, M2E21

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

        oA00_A = oA00.assembled
        oA00_A = oA00_A.___PRIVATE_gather_M_to_core___(clean_local=True)
        oA00_A = GlobalMatrix(oA00_A)

        oB_0_A = oB_0.assembled
        oB_0_A = oB_0_A.___PRIVATE_gather_V_to_core___(clean_local=True)
        oB_0_A = GlobalVector(oB_0_A)

        if rAnk == mAster_rank:

            M0_ = oA._M_[0:u2.GLOBAL_num_dofs]
            oA._M_ = oA._M_[u2.GLOBAL_num_dofs:]
            M01 = M0_[:, u2.GLOBAL_num_dofs:]
            M0_ = spspa.hstack((oA00_A.M, M01), format='csr')
            oA._M_ = spspa.vstack((M0_, oA._M_), format='csc')

            # warnings.filterwarnings("ignore")
            # oA._M_[0:u2.GLOBAL_num_dofs, 0:u2.GLOBAL_num_dofs] = oA00_A.M
            # warnings.filterwarnings("default")

            ob.V[0:u2.GLOBAL_num_dofs] = oB_0_A.V
        del oB_0_A, oA00_A

        if tk == t0:  # first step
            X0_0 = u2.cochain.globe
            X0_1 = w1.cochain.globe
            X0_2 = DistributedVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
            X0 = TLO.concatenate((X0_0, X0_1, X0_2))
        else:
            X0 = OUT_R[0]
        oR, _, _, _, mo = getattr(scipy_sparse_linalg, solver)(
                        oA, ob, X0, tol=tol, restart=restart, maxiter=maxiter)
        OUT_R[0] = oR
        oR.DO_distribute_to(u2, w1, P3)

        du2 = u2.coboundary()
        du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
        du2.TW.current_time = tk1
        du2.TW.___DO_push_all_to_instant___()
        DIV_L2_error_tk1 = du2.error.L()
        KE2_tk1 = 0.5 * u2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume
        E1_tk1 = 0.5 * w1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume

        # ... inner
        iA00_A = iA00.assembled
        iA00_A = iA00_A.___PRIVATE_gather_M_to_core___(clean_local=True)
        iA00_A = GlobalMatrix(iA00_A)

        iB_0_A = iB_0.assembled
        iB_0_A = iB_0_A.___PRIVATE_gather_V_to_core___(clean_local=True)
        iB_0_A = GlobalVector(iB_0_A)

        if rAnk == mAster_rank:
            ____ = iA._M_[0:u1.GLOBAL_num_dofs]
            iA._M_ = iA._M_[u1.GLOBAL_num_dofs:]
            ____ = ____[:, u1.GLOBAL_num_dofs:]
            ____ = spspa.hstack((iA00_A.M, ____), format='csr')
            iA._M_ = spspa.vstack((____, iA._M_), format='csr')

            # warnings.filterwarnings("ignore")
            # iA._M_[0:u1.GLOBAL_num_dofs, 0:u1.GLOBAL_num_dofs] = iA00_A.M
            # warnings.filterwarnings("default")

            ib.V[0:u1.GLOBAL_num_dofs] = iB_0_A.V
        del iB_0_A, iA00_A

        X0 = INN_R[0]
        iR, _, _, _, mi = getattr(scipy_sparse_linalg, solver)(
                        iA, ib, X0, tol=tol, restart=restart, maxiter=maxiter)
        INN_R[0] = iR
        _u1_old_cochain_ = u1.cochain.local
        iR.DO_distribute_to(u1, P0)

        _u1_new_cochain_ = u1.cochain.local

        mean_u1_cochain_local_at_tk = dict()
        for i in _u1_old_cochain_:
            mean_u1_cochain_local_at_tk[i] = (_u1_old_cochain_[i] + _u1_new_cochain_[i]) / 2

        u1.cochain.local = mean_u1_cochain_local_at_tk  # we then have u1 cochain @ tk
        KE1_tk1 = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume
        H1_tk1 = u1.DO.compute_L2_inner_product_energy_with(w1, M=M1)
        u1u2_diff_tk1 = u2.DO.compute_L2_diff_from(u1)

        w2.cochain.local = u1.coboundary.cochain_local
        H2_tk1 = u2.DO.compute_L2_inner_product_energy_with(w2, M=M2)
        E2_tk1 = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume
        w1w2_diff_tk1 = w2.DO.compute_L2_diff_from(w1)

        if save_uw:
            if check_multiple_close(tk1, 0.1) and check_almost_in_range(tk1, 8.7, 9.5):
                TK1 = round(tk1, 1)
                save([u1, u2, w1, w2], f'UUWW_TGV_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            elif check_multiple_close(tk1, 1):
                TK1 = round(tk1)
                save([u1, u2, w1, w2], f'UUWW_TGV_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            else:
                pass
        u1.cochain.local = _u1_new_cochain_  # renew u1 cochain to time tk+half
        w2.cochain.local = u1.coboundary.cochain_local  # renew w2 cochain to time tk+half
        KE1_tk1h = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1) / Volume
        E2_tk1h = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2) / Volume
        if save_uw:
            if check_multiple_close(tk1, 0.1) and check_almost_in_range(tk1, 8.3, 9.5):
                TK1 = round(tk1, 1)
                save([u1, w2], f'UWih_TGV_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            elif check_multiple_close(tk1, 1):
                TK1 = round(tk1)
                save([u1, w2], f'UWih_TGV_Re{Re}_N{N}k{k}t{t}Steps{steps}_t{TK1}')
            else:
                pass

        message = [f'ITERATION cost: {int((time() - ts) * 100) / 100}',
                   'Inner solver: ' + mi,
                   'Outer solver: ' + mo]

        # print(KE1_tk1, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1)

        return 1, 0, message, KE1_tk1, KE1_tk1h, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, E2_tk1h, \
               u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1

    # SOLVER(t0, t0+dt)

    SI(SOLVER, [KE1_t0, KE1_t0h, KE2_t0, H1_t0, H2_t0, E1_t0, E2_t0, E2_t0h, u1u2_diff_t0,
                w1w2_diff_t0, DIV_L2_error_t0])
    SI.run()

    return SI

if __name__ == '__main__':
    # mpiexec python _3dCSCG\TESTS\__unittest_scripts__\icpsNS_TGV_no_hybrid.py

    NoHy_TGV(N=2, k=3, t=1, steps=20, Re=500, tol=1e-3, restart=30, maxiter=10, solver='gcrotmk', save_uw=False)