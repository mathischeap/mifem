

import sys
if './' not in sys.path: sys.path.append('../__unittest_scripts__/')

from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.data_structures.global_matrix.main import GlobalMatrix, GlobalVector
from scipy import sparse as spspa
from tools.iterators.simple import SimpleIterator
import tools.linear_algebra.deprecated.operators as TLO
import tools.linear_algebra.solvers.serial.deprecated as scipy_sparse_linalg
from time import time

from root.config.main import *

def manu_conserving_solver(N, k, t, steps):
       # ... define the problem ...
       t0 = 0
       dt = (t-t0)/steps
       # nu = 0
       # rho = 1
       quad_degree = [2*N, 2*N, 2*N]
       RDF_filename = f'RDF_manuC_N{N}K{k}t{t}Steps{steps}'
       name         = f'ITR_manuC_N{N}K{k}t{t}Steps{steps}'
       auto_save_frequency = 5
       monitor_factor = 1

       SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                           auto_save_frequency=auto_save_frequency,
                           monitor_factor=monitor_factor,
                           RDF_filename=RDF_filename,
                           name=name)

       mesh = MeshGenerator('crazy_periodic', c=0.0)([k, k, k])
       space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N), ('Lobatto', N)])
       FC = FormCaller(mesh, space)
       es = ExactSolutionSelector(mesh)('icpsNS:sincosRC')
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
       CP1 = w1.special.cross_product_1f__ip_1f(u1, u1, quad_degree=quad_degree)
       CP2 = w2.special.cross_product_2f__ip_2f(u2, u2, quad_degree=quad_degree)

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
       KE1_t0 = 0.5 * u1.do.compute_L2_energy_with(M=M1)
       KE2_t0 = 0.5 * u2.do.compute_L2_energy_with(M=M2)
       H1_t0  =       u1.do.compute_L2_energy_with(w1, M=M1)
       H2_t0  =       u2.do.compute_L2_energy_with(w2, M=M2)
       E1_t0  = 0.5 * w1.do.compute_L2_energy_with(M=M1)
       E2_t0  = 0.5 * w2.do.compute_L2_energy_with(M=M2)
       du2 = u2.coboundary()
       du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
       du2.TW.current_time = t0
       du2.TW.___DO_push_all_to_instant___()
       DIV_L2_error_t0 = du2.error.L()
       DIV_L_inf_error_t0 = du2.error.L(n='infinity')
       u1u2_diff_t0 = u2.do.compute_Ln_diff_from(u1)
       w1w2_diff_t0 = w2.do.compute_Ln_diff_from(w1)

       # set up inner half integer time step systems ......
       lhs00_Hf = 2*M1/dt + 0.5*CP1
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
       iA = TLO.bmat(lhs)
       iA = iA.___PRIVATE_gather_M_to_core___(clean_local=True)
       iA = GlobalMatrix(iA)
       assert iA.IS.master_dominating
       if rAnk != mAster_rank: iA._M_ = None

       B0 = (2 * M1 / dt - 0.5 * CP1) @ u1.cochain.EWC
       B0.gathering_matrix = u1
       B0 = B0.assembled
       B1 = GlobalVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
       ib = TLO.concatenate([B0, B1])
       ib = ib.___PRIVATE_gather_V_to_core___(clean_local=True)
       ib = GlobalVector(ib)
       assert ib.IS.master_dominating

       del lhs00_Hf, lhs00_Hf_A, M1E10, M1E10_A, E01M1_A, lhs11_A, lhs, B0, B1

       iR = scipy_sparse_linalg.spsolve(iA, ib)[0]
       iR.___PRIVATE_be_distributed_to___(u1, P0)
       w2.cochain.local = u1.coboundary.cochain_local
       KE1_t0h = 0.5 * u1.do.compute_L2_energy_with(M=M1)
       E2_t0h  = 0.5 * w2.do.compute_L2_energy_with(M=M2)

       u1_p_f2 = u1.special.___PRIVATE_projected_into_2form_exactly___()
       D_u1_p_f2 = u1_p_f2.coboundary()
       D_u1_p_f2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
       D_u1_p_f2.TW.current_time = dt / 2
       D_u1_p_f2.TW.___DO_push_all_to_instant___()
       DIV_u1_L2_error_tk0h = D_u1_p_f2.error.L()
       DIV_u1_L_inf_error_tk0h = D_u1_p_f2.error.L(n='infinity', quad_density=1000000)



       iA00 =  M1/dt + 0.5*CP1
       iA00.gathering_matrices = (u1, u1)
       iB_0 = (M1/dt - 0.5*CP1) @ u1.cochain.EWC
       iB_0.gathering_matrix = u1

       # set up outer integer time step systems ......
       oA00 = M2 / dt + 0.5 * CP2
       oA00.gathering_matrices = (u2, u2)

       E23M3 = E23 @ M3
       E23M3.gathering_matrices = (u2, P3)
       mE23M3_A = - E23M3.assembled

       E12M2 = E12 @ M2
       E12M2.gathering_matrices = (w1, u2)
       mE12M2_A = - E12M2.assembled

       M1.gathering_matrices = (w1, w1)
       M1_A = M1.assembled

       E32.gathering_matrices = (P3, u2)
       E32_A = E32.assembled


       lhs = [[ None   , None, mE23M3_A],  # u2
              [mE12M2_A, M1_A,  None   ],  # w1
              [ E32_A  , None,  None   ]]

       oA = TLO.bmat(lhs)
       oA = oA.___PRIVATE_gather_M_to_core___(clean_local=True)
       oA = GlobalMatrix(oA)
       assert oA.IS.master_dominating
       if rAnk != mAster_rank: oA._M_ = None

       oB_0 = (M2 / dt - 0.5 * CP2) @ u2.cochain.EWC
       oB_0.gathering_matrix = u2
       B0 = oB_0.assembled
       B1 = GlobalVector(spspa.csc_matrix((w1.GLOBAL_num_dofs, 1)))
       B2 = GlobalVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
       ob = TLO.concatenate([B0, B1, B2])
       ob = ob.___PRIVATE_gather_V_to_core___(clean_local=True)
       ob = GlobalVector(ob)
       assert ob.IS.master_dominating
       del E23M3, mE23M3_A, E12M2, mE12M2_A, M1_A, E32_A, lhs, B0, B1, B2


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
              DIV_L_inf_error :
              DIV_u1_L2_error_tk1h :
              DIV_u1_L_inf_error_tk1h :
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
                  ____ = oA._M_[0:u2.GLOBAL_num_dofs]
                  oA._M_ = oA._M_[u2.GLOBAL_num_dofs:]
                  ____ = ____[:, u2.GLOBAL_num_dofs:]
                  ____ = spspa.hstack((oA00_A.M, ____), format='csr')
                  oA._M_ = spspa.vstack((____, oA._M_), format='csc')
                  ob.V[0:u2.GLOBAL_num_dofs] = oB_0_A.V
              del oB_0_A, oA00_A

              oR, _, _, _, mo = scipy_sparse_linalg.spsolve(oA, ob)
              oR.___PRIVATE_be_distributed_to___(u2, w1, P3)

              du2 = u2.coboundary()
              du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
              du2.TW.current_time = tk1
              du2.TW.___DO_push_all_to_instant___()
              DIV_L2_error_tk1 = du2.error.L()
              DIV_L_inf_error_tk1 = du2.error.L(n = 'infinity')
              KE2_tk1 = 0.5 * u2.do.compute_L2_energy_with(M=M2)
              E1_tk1 = 0.5 * w1.do.compute_L2_energy_with(M=M1)

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
                  ib.V[0:u1.GLOBAL_num_dofs] = iB_0_A.V
              del iB_0_A, iA00_A

              iR, _, _, _, mi = scipy_sparse_linalg.spsolve(iA, ib)
              _u1_old_cochain_ = u1.cochain.local
              iR.___PRIVATE_be_distributed_to___(u1, P0)

              _u1_new_cochain_ = u1.cochain.local

              mean_u1_cochain_local_at_tk = dict()
              for i in _u1_old_cochain_:
                  mean_u1_cochain_local_at_tk[i] = (_u1_old_cochain_[i] + _u1_new_cochain_[i]) / 2

              u1.cochain.local = mean_u1_cochain_local_at_tk  # we then have u1 cochain @ tk
              KE1_tk1 = 0.5 * u1.do.compute_L2_energy_with(M=M1)
              H1_tk1 = u1.do.compute_L2_energy_with(w1, M=M1)
              u1u2_diff_tk1 = u2.do.compute_Ln_diff_from(u1)

              w2.cochain.local = u1.coboundary.cochain_local
              H2_tk1 = u2.do.compute_L2_energy_with(w2, M=M2)
              E2_tk1 = 0.5 * w2.do.compute_L2_energy_with(M=M2)
              w1w2_diff_tk1 = w2.do.compute_Ln_diff_from(w1)

              u1.cochain.local = _u1_new_cochain_  # renew u1 cochain to time tk+half
              w2.cochain.local = u1.coboundary.cochain_local  # renew w2 cochain to time tk+half
              KE1_tk1h = 0.5 * u1.do.compute_L2_energy_with(M=M1)
              E2_tk1h = 0.5 * w2.do.compute_L2_energy_with(M=M2)

              u1_p_f2 = u1.special.___PRIVATE_projected_into_2form_exactly___()
              D_u1_p_f2 = u1_p_f2.coboundary()
              D_u1_p_f2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
              D_u1_p_f2.TW.current_time = tk1 + dt / 2
              D_u1_p_f2.TW.___DO_push_all_to_instant___()
              DIV_u1_L2_error_tk1h = D_u1_p_f2.error.L()
              DIV_u1_L_inf_error_tk1h = D_u1_p_f2.error.L(n='infinity', quad_density=1000000)


              message = [f'ITERATION cost: {int((time()-ts)*100)/100}',
                         'Inner solver: ' + mi,
                         'Outer solver: ' + mo]



              # print(KE1_tk1, KE1_tk1h, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, E2_tk1h,
              #       u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1)

              return 1, 0, message, KE1_tk1, KE1_tk1h, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, E2_tk1h, \
                     u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1, DIV_L_inf_error_tk1, \
                     DIV_u1_L2_error_tk1h, DIV_u1_L_inf_error_tk1h


       SI(SOLVER, [KE1_t0, KE1_t0h, KE2_t0, H1_t0, H2_t0, E1_t0, E2_t0, E2_t0h, u1u2_diff_t0,
                   w1w2_diff_t0, DIV_L2_error_t0, DIV_L_inf_error_t0,
                   DIV_u1_L2_error_tk0h, DIV_u1_L_inf_error_tk0h
                   ])
       SI.run()

       return SI


if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG/TESTS/__unittest_scripts__/icpsNS_manu_conserving_solver.py

    manu_conserving_solver(2, 3, 10, 200)