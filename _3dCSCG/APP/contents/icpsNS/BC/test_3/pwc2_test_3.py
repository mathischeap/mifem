# -*- coding: utf-8 -*-
"""
$ mpiexec python _3dCSCG\APP\contents\icpsNS\BC\test_3\pwc2_test_3.py

"""


import sys

if './' not in sys.path: sys.path.append('./')

from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from TOOLS.__DEPRECATED__.linear_system.main import LinearSystem
from TOOLS.linear_algebra.data_structures import GlobalMatrix, GlobalVector, DistributedVector
from TOOLS.iterator import SimpleIterator
from scipy import sparse as spspa
import TOOLS.linear_algebra.__DEPRECATED__.operators as TLO
from root.mifem import save, read




t0 = 0
t = 2
steps=200
start_with = None
half_step_u1w2 = None



# do not change below ....
N = 3
element_layout = [25, 15, [1,2,4,8,16,32,64,64,64,64,64,64,64,64,64,64,32,16,8,4,2,1]]
F = 1
nu = 0.01
rho = 1
restart = 300
tol = 1e-5
maxiter = 50

t0 = int(t0)
t = int(t)

conserving = False
RDF_filename = f'RDF_PWC2_TEST_3_t0{int(t0)}_t{int(t)}'
name         = f'Iterator_PWC2_TEST_3_t0{int(t0)}_t{int(t)}'
auto_save_frequency = 5
monitor_factor = 1

dt = (t-t0) / steps

SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                    auto_save_frequency=auto_save_frequency,
                    monitor_factor=monitor_factor,
                    RDF_filename=RDF_filename,
                    name=name)

quad_degree = [2*N, 2*N, 2*N]


mesh = MeshGenerator('pwc', l=1, w=0.5, h=0.5)(element_layout)

space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N), ('Lobatto', N)])
FC = FormCaller(mesh, space)
es = ExactSolutionSelector(mesh)('icpsNS:CBFx', f=F, nu=nu, rho=rho)

P0 = FC('0-f', is_hybrid=False, orientation='inner', name='inner-total-pressure')
u1 = FC('1-f', is_hybrid=False, orientation='inner', name='inner-velocity')
w2 = FC('2-f', is_hybrid=False, orientation='inner', name='inner-vorticity')
w1 = FC('1-f', is_hybrid=False, orientation='outer', name='outer-vorticity')
u2 = FC('2-f', is_hybrid=True, orientation='outer', name='outer-velocity')
P3 = FC('3-f', is_hybrid=True, orientation='outer', name='outer-total-pressure')
t2 = FC('2-t', orientation='outer', name='Trace2_outer')

f = es.status.body_force
f.current_time = t0
FI = f.DO.inner_product_with_space_of(u1, quad_degree=[N + 1, N + 1, N + 1])
FO = f.DO.inner_product_with_space_of(u2, quad_degree=[N + 1, N + 1, N + 1])

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
N2 = t2.coboundary.trace_matrix

if start_with is None:
    u1.TW.func.DO.set_func_body_as(es.status.velocity)
    u2.TW.func.DO.set_func_body_as(es.status.velocity)
    w1.TW.func.DO.set_func_body_as(es.status.vorticity)
    w2.TW.func.DO.set_func_body_as(es.status.vorticity)
    u1.TW.current_time = t0
    w2.TW.current_time = t0
    w1.TW.current_time = t0
    u2.TW.current_time = t0
    u1.TW.___DO_push_all_to_instant___()
    w2.TW.___DO_push_all_to_instant___()
    w1.TW.___DO_push_all_to_instant___()
    u2.TW.___DO_push_all_to_instant___()
    u1.discretize()
    w2.discretize()
    w1.discretize()
    u2.discretize()
else:
    U1, U2, W1, W2 = read(start_with)
    u1.DO.resemble(U1)
    u2.DO.resemble(U2)
    w1.DO.resemble(W1)
    w2.DO.resemble(W2)

KE1_t0 = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1)
KE2_t0 = 0.5 * u2.DO.compute_L2_inner_product_energy_with(M=M2)
H1_t0  =       u1.DO.compute_L2_inner_product_energy_with(w1, M=M1)
H2_t0  =       u2.DO.compute_L2_inner_product_energy_with(w2, M=M2)
E1_t0  = 0.5 * w1.DO.compute_L2_inner_product_energy_with(M=M1)
E2_t0  = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2)

du2 = u2.coboundary()
du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
du2.TW.current_time = t0
du2.TW.___DO_push_all_to_instant___()
DIV_L2_error_t0 = du2.error.L()
u1u2_diff_t0 = u2.DO.compute_L2_diff_from(u1)
w1w2_diff_t0 = w2.DO.compute_L2_diff_from(w1)

E12M2E21 = E12 @ M2 @ E21
lhs00_Hf = 2*M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21

lhs00_Hf.gathering_matrices = (u1, u1)
lhs00_Hf_A = lhs00_Hf.assembled
del lhs00_Hf

M1E10 = M1 @ E10
M1E10.gathering_matrices = (u1, P0)

M1E10_A = M1E10.assembled
E01M1_A = M1E10_A.T
E01M1_A.DO.clear_row(-1)
lhs11_A = GlobalMatrix((P0.GLOBAL_num_dofs, P0.GLOBAL_num_dofs))
lhs11_A.DO.set_value(-1, -1, 1)
del M1E10

# no-slip BC for u1 ...
u1_WALL_DoFs = u1.numbering.GLOBAL_boundary_dofs_ravel
M1E10_A.DO.clear_rows(u1_WALL_DoFs)
lhs00_Hf_A.DO.identify_rows(u1_WALL_DoFs)

lhs = [[ lhs00_Hf_A, M1E10_A],
       [-E01M1_A   , lhs11_A]]
del E01M1_A, M1E10_A, lhs11_A, lhs00_Hf_A

B0 = (2*M1/dt - 0.5*CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC + FI
B0.gathering_matrix = u1
B0 = B0.assembled

B0.DO.clear_entries(u1_WALL_DoFs)

B1 = GlobalVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
LS_hal = LinearSystem(lhs, rhs=[B0, B1])

x0_0 = u1.cochain.globe
x0_1 = DistributedVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
x0 = TLO.concatenate((x0_0, x0_1))

if half_step_u1w2 is None:
    LS_hal.solve('GMRES', 'MPI1')(x0, restart=restart, tol=tol, maxiter=maxiter, conserving=conserving)
    LS_hal.solve.results.DO_distribute_to(u1, P0)
    w2.cochain.local = u1.coboundary.cochain_local
else:
    U1, W2 = read(half_step_u1w2)
    u1.DO.resemble(U1)
    w2.DO.resemble(W2)

# KE1_t0h = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1)
# E2_t0h  = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2)

LS_half_lhs00 =  M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21                         # assemble it before using
LS_half_B0    = (M1/dt - 0.5*CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC  + FI # assemble it before using
LS_half_lhs00.gathering_matrices = (u1, u1)
LS_half_B0.gathering_matrix = u1

lhs00 = M2/dt + 0.5*CP2
E23M3 = E23 @ M3
M2E21 = M2 @ E21
M2E21.gathering_matrices = (u2, w1)
M2E21_A = M2E21.assembled
N2.gathering_matrices = (t2, u2)
N2_A = N2.assembled
M1.gathering_matrices=(w1, w1)
M1_A = M1.assembled
E32.gathering_matrices = (P3, u2)


t2_WALL_DoFs = t2.numbering.GLOBAL_boundary_dofs_ravel
for i in range(t2.GLOBAL_num_dofs):
    if i not in t2_WALL_DoFs:
        break
lhs33_A = GlobalMatrix((t2.GLOBAL_num_dofs, t2.GLOBAL_num_dofs))
lhs33_A.DO.set_value(i, i, 1)
N2_A.DO.clear_row(i)


lhs = [[ lhs00    ,-E23M3, 0.5*nu*M2E21_A, N2_A.T ], # u2
       [ E32      , None , None          , None   ], # P3
       [-M2E21_A.T, None , M1_A          , None   ], # w1
       [ N2_A     , None , None          , lhs33_A]] # t2
del M2E21_A, E23M3, lhs00

# E23M3.gathering_matrices = (u2, P3)
# E23M3_A = E23M3.assembled
# E32_A = E32.assembled
# del E23M3
# lhs00.gathering_matrices = (u2, u2)
# lhs = [[ lhs00    ,-E23M3_A, 0.5*nu*M2E21_A, N2_A.T ], # u2
#        [ E32_A    , None   , None          , None   ], # P3
#        [-M2E21_A.T, None   , M1_A          , None   ], # w1
#        [ N2_A     , None   , None          , lhs33_A]] # t2


B0 = (M2/dt - 0.5*CP2) @ u2.cochain.EWC - 0.5*nu*M2E21 @ w1.cochain.EWC + FO
B0.gathering_matrix = u2
B1 = GlobalVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
B2 = GlobalVector(spspa.csc_matrix((w1.GLOBAL_num_dofs, 1)))
B3 = GlobalVector(spspa.csc_matrix((t2.GLOBAL_num_dofs, 1)))
LS_int = LinearSystem(lhs, rhs=[B0, B1, B2, B3])

del M2E21, N2_A, M1_A, lhs33_A, E12M2E21, E32

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
    KE2_tk :
    H1_tk :
    H2_tk :
    E1_tk :
    E2_tk :
    u2u1_L2_diff :
    w2w1_L2_diff :
    DIV_L2_error :
    """
    assert tk1 == tk + dt


    if tk == t0: # first step
        X0_0 = w1.cochain.globe
        X0_1 = DistributedVector(spspa.csc_matrix((t2.GLOBAL_num_dofs, 1)))
        X0 = TLO.concatenate((X0_0, X0_1))
    else:
        X0 = LS_int.solve.others[0]
    LS_int.solve('icpsNS', 'OS1')(X0, restart=restart, tol=tol, maxiter=maxiter, conserving=conserving)

    # if tk == t0: # first step
    #     X0_0 = u2.cochain.globe
    #     X0_1 = DistributedVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
    #     X0_2 = w1.cochain.globe
    #     X0_3 = DistributedVector(spspa.csc_matrix((t2.GLOBAL_num_dofs, 1)))
    #     X0 = TLO.concatenate((X0_0, X0_1, X0_2, X0_3))
    # else:
    #     X0 = LS_int.solve.results
    # LS_int.lhs[0][0] = lhs00.assembled
    # LS_int.rhs[0]    = B0.assembled
    # LS_int.solve('GMRES', 'MPI1')(X0, restart=restart, tol=tol, maxiter=maxiter, conserving=conserving)


    LS_int.solve.results.DO_distribute_to(u2, P3, w1, t2)

    du2 = u2.coboundary()
    du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
    du2.TW.current_time = tk1
    du2.TW.___DO_push_all_to_instant___()
    DIV_L2_error_tk1 = du2.error.L()
    KE2_tk1 = 0.5 * u2.DO.compute_L2_inner_product_energy_with(M=M2)
    E1_tk1  = 0.5 * w1.DO.compute_L2_inner_product_energy_with(M=M1)

    # ...
    LS_hal_LHS00 = LS_half_lhs00.assembled
    LS_hal_RHS00 = LS_half_B0.assembled
    LS_hal_LHS00.DO.identify_rows(u1_WALL_DoFs)
    LS_hal_RHS00.DO.clear_entries(u1_WALL_DoFs)
    LS_hal.lhs[0][0] = LS_hal_LHS00
    LS_hal.rhs[0] = LS_hal_RHS00

    if tk == t0:
        if half_step_u1w2 is None:
            X0 = LS_hal.solve.results
        else:
            x0_0 = u1.cochain.globe
            x0_1 = DistributedVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
            X0 = TLO.concatenate((x0_0, x0_1))
    else:
        X0 = LS_hal.solve.results

    LS_hal.solve('GMRES', 'MPI1')(X0, restart=restart, tol=tol, maxiter=maxiter, conserving=conserving)
    _u1_old_cochain_ = u1.cochain.local
    LS_hal.solve.results.DO_distribute_to(u1)
    _u1_new_cochain_ = u1.cochain.local

    mean_u1_cochain_local_at_tk = dict()
    for i in _u1_old_cochain_:
        mean_u1_cochain_local_at_tk[i] = (_u1_old_cochain_[i] + _u1_new_cochain_[i]) / 2

    u1.cochain.local = mean_u1_cochain_local_at_tk # we then have u1 cochain @ tk
    KE1_tk1 = 0.5 * u1.DO.compute_L2_inner_product_energy_with(M=M1)
    H1_tk1  =       u1.DO.compute_L2_inner_product_energy_with(w1, M=M1)
    u1u2_diff_tk1 = u2.DO.compute_L2_diff_from(u1)

    w2.cochain.local = u1.coboundary.cochain_local
    H2_tk1  =       u2.DO.compute_L2_inner_product_energy_with(w2, M=M2)
    E2_tk1  = 0.5 * w2.DO.compute_L2_inner_product_energy_with(M=M2)
    w1w2_diff_tk1 = w2.DO.compute_L2_diff_from(w1)

    if abs(round(tk1) - tk1) < 1e-8: save([u1, u2, w1, w2],
        f'PWC2_uuww_TEST_3_t{int(round(tk1))}')

    u1.cochain.local = _u1_new_cochain_              # renew u1 cochain to time tk+half
    w2.cochain.local = u1.coboundary.cochain_local   # renew w2 cochain to time tk+half

    if abs(round(tk1) - tk1) < 1e-8: save([u1, w2],
        f'PWC2_u1w2_TEST_3_t{int(round(tk1))}h')

    message = ['Inner solver: ' + LS_hal.solve.message, 'Outer solver: ' + LS_int.solve.message]

    return 1, 0, message, KE1_tk1, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, \
           u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1


SI(SOLVER, [KE1_t0, KE2_t0, H1_t0, H2_t0, E1_t0, E2_t0, u1u2_diff_t0,
            w1w2_diff_t0, DIV_L2_error_t0])
SI.run()