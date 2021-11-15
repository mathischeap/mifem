# -*- coding: utf-8 -*-
"""
$ mpiexec python _3dCSCG\APP\contents\icpsNS\BC\pwc1.py

"""
import sys
if './' not in sys.path: sys.path.append('./')
# from config import *
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from TOOLS.__DEPRECATED__.linear_system.main import LinearSystem
from TOOLS.linear_algebra.data_structures import GlobalMatrix, GlobalVector, DistributedVector
from root.mifem import read, save
import TOOLS.linear_algebra.__DEPRECATED__.operators as TLO
from scipy import sparse as spspa
from TOOLS.iterator import SimpleIterator

N = 3
k = 8
t0 = 0
t = 1
steps = 400
# hsP0 = 'pwc1_hsP0_N1k2elements20.mi'
hsP0 = None
hsP3 = None
save_HS = False

restart = 150
maxiter =50
tol = 1e-6

F = 1
nu = 0.1
rho = 1

dt = (t-t0)/steps
quad_degree = [2 * N, 2 * N, 2 * N]

auto_save_frequency = 5
monitor_factor = 1
RDF_filename = 'RDF_PWC1_t{}N{}K{}Steps{}'.format(t,N, k, steps)
name = 'Iterator_PWC1_t{}N{}K{}Steps{}'.format(t, N, k, steps)

SI = SimpleIterator(t0=t0, dt=dt, max_steps=steps,
                    auto_save_frequency=auto_save_frequency,
                    monitor_factor=monitor_factor,
                    RDF_filename=RDF_filename,
                    name=name)

# mesh = MeshGenerator('pwc', l=1, w=1, h=1)([k, k, [1,2,4,8,16,16,16,8,4,2,1]])
mesh = MeshGenerator('pwc', l=1, w=1, h=1)([2, 2, [1,2,1]])
space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N), ('Lobatto', N)])
FC = FormCaller(mesh, space)
es = ExactSolutionSelector(mesh)('icpsNS:CBFx', f=F, nu=nu, rho=rho)

P0 = FC('0-f', is_hybrid=False, orientation='inner', name='inner-total-pressure')
u1 = FC('1-f', is_hybrid=False, orientation='inner', name='inner-velocity')
w2 = FC('2-f', is_hybrid=False, orientation='inner', name='inner-vorticity')
w1 = FC('1-f', is_hybrid=False, orientation='outer', name='outer-vorticity')
u2 = FC('2-f', is_hybrid=False, orientation='outer', name='outer-velocity')
P3 = FC('3-f', is_hybrid=False, orientation='outer', name='outer-total-pressure')

f = es.status.body_force

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


# ...
u1.TW.func.___DO_set_func_body_as___(es.status.velocity)
u2.TW.func.___DO_set_func_body_as___(es.status.velocity)
w1.TW.func.___DO_set_func_body_as___(es.status.vorticity)
w2.TW.func.___DO_set_func_body_as___(es.status.vorticity)

u1.TW.current_time = t0
w2.TW.current_time = t0
u1.TW.___DO_push_all_to_instant___()
w2.TW.___DO_push_all_to_instant___()
u1.discretize()
w2.discretize()

w1.TW.current_time = t0
u2.TW.current_time = t0
w1.TW.___DO_push_all_to_instant___()
u2.TW.___DO_push_all_to_instant___()
w1.discretize()
u2.discretize()


# ...
E12M2E21 = E12 @ M2 @ E21
lhs00_Hf = 2 * M1 / dt + 0.5 * CP1 + 0.5 * nu * E12M2E21

lhs00_Hf.gathering_matrices = (u1, u1)
lhs00_Hf_A = lhs00_Hf.assembled
del lhs00_Hf

M1E10 = M1 @ E10
M1E10.gathering_matrices = (u1, P0)
M1E10_A = M1E10.assembled
E01M1_A = M1E10_A.T

# fix the pressure at one dof ...
E01M1_A.DO.clear_row(0)
lhs11_A = GlobalMatrix((P0.GLOBAL_num_dofs, P0.GLOBAL_num_dofs))
lhs11_A.DO.set_value(0, 0, 1)

# no-slip BC for u1 ...
u1_WALL_DoFs = u1.numbering.GLOBAL_boundary_dofs_ravel
M1E10_A.DO.clear_rows(u1_WALL_DoFs)
lhs00_Hf_A.DO.identify_rows(u1_WALL_DoFs)

# lhs ...
lhs = [[ lhs00_Hf_A, M1E10_A],
       [-E01M1_A   , lhs11_A]]

FI = f.___PRIVATE_do_inner_product_with_space_of___(u1, quad_degree=[N + 1, N + 1, N + 1])
FO = f.___PRIVATE_do_inner_product_with_space_of___(u2, quad_degree=[N + 1, N + 1, N + 1])

f.current_time = t0 + dt / 4
B0 = (2 * M1 / dt - 0.5 * CP1 - 0.5 * nu * E12M2E21) @ u1.cochain.EWC + FI
B0.gathering_matrix = u1
B0 = B0.assembled
B0.DO.clear_entries(u1_WALL_DoFs)
B1 = GlobalVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
LS_hal = LinearSystem(lhs, rhs=[B0, B1])


x0_0 = u1.cochain.globe
if hsP0 is None:
    x0_1 = DistributedVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
else:
    temp_P0 = read(hsP0)
    P0.DO.resemble(temp_P0)
    x0_1 = P0.cochain.globe

x0 = TLO.concatenate((x0_0, x0_1))
LS_hal.solve('GMRES', 'MPI1')(x0, conserving=False,
        restart=restart, tol=tol, maxiter=maxiter)

LS_hal.solve.results.DO_distribute_to(u1, P0)
P0.TW.current_time = t0 + dt/4
if save_HS: save(P0, f'pwc1_hsP0_N{N}k{k}elements{mesh.elements.GLOBAL_num}')

w2.cochain.local = u1.coboundary.cochain_local
LS_half_lhs00 =  M1/dt + 0.5*CP1 + 0.5*nu*E12M2E21                         # assemble it before using
LS_half_B0    = (M1/dt - 0.5*CP1 - 0.5*nu*E12M2E21) @ u1.cochain.EWC + FI  # assemble it before using
LS_half_lhs00.gathering_matrices = (u1, u1)
LS_half_B0.gathering_matrix = u1





lhs00_Is = M2 / dt + 0.5 * CP2
lhs00_Is.gathering_matrices = (u2, u2)

M3E32 = M3 @ E32
M3E32.gathering_matrices = (P3, u2)
M3E32_A = M3E32.assembled
E23M3_A = M3E32_A.T


M2E21 = M2 @ E21
M2E21.gathering_matrices = (u2, w1)
M2E21_A = M2E21.assembled

E12M2 = E12 @ M2
E12M2.gathering_matrices = (w1, u2)
E12M2_A = E12M2.assembled

M1.gathering_matrices = (w1, w1)
M1_A = M1.assembled

# wall BC for u2 ...
u2_WALL_DoFs = u2.numbering.GLOBAL_boundary_dofs_ravel
M2E21_A.DO.clear_rows(u2_WALL_DoFs)
E23M3_A.DO.clear_rows(u2_WALL_DoFs)

# fix the pressure at one dof ...
lhs22_A = GlobalMatrix((P3.GLOBAL_num_dofs, P3.GLOBAL_num_dofs))
M3E32_A.DO.clear_row(0)
lhs22_A.DO.set_value(0, 0, 1)


lhs = [[ lhs00_Is, 0.5 * nu * M2E21_A, -E23M3_A],  # u2
       [-E12M2_A ,  M1_A             ,     None],  # w1
       [ M3E32_A ,  None             ,  lhs22_A]]  # P3

LS_int_B0 = (M2/dt - 0.5*CP2) @ u2.cochain.EWC - 0.5*nu*M2E21 @ w1.cochain.EWC + FO
LS_int_B0.gathering_matrix = u2
B1 = GlobalVector(spspa.csc_matrix((w1.GLOBAL_num_dofs, 1)))
B2 = GlobalVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))

LS_int = LinearSystem(lhs, rhs=[LS_int_B0, B1, B2])

del M2E21


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
        X0_0 = u2.cochain.globe
        X0_1 = w1.cochain.globe

        if hsP3 is None:
            X0_2 = DistributedVector(spspa.csc_matrix((P3.GLOBAL_num_dofs, 1)))
        else:
            temp_P3 = read(hsP3)
            P3.DO.resemble(temp_P3)
            X0_2 = P3.cochain.globe

        X0 = TLO.concatenate((X0_0, X0_1, X0_2))
    else:
        X0 = LS_int.solve.results

    LHS_INT_00_A = lhs00_Is.assembled
    RHS_INT_00_A = LS_int_B0.assembled
    LHS_INT_00_A.DO.identify_rows(u2_WALL_DoFs)
    RHS_INT_00_A.DO.clear_entries(u2_WALL_DoFs)
    LS_int.lhs[0][0] = LHS_INT_00_A
    LS_int.rhs[0]    = RHS_INT_00_A

    print(LS_int.condition.condition_number, flush=True)

    LS_int.solve('GMRES', 'MPI1')(X0, restart=restart, tol=tol, maxiter=maxiter, conserving=False)

    LS_int.solve.results.DO_distribute_to(u2, w1, P3)
    du2 = u2.coboundary()
    du2.TW.func.___DO_set_func_body_as___(es.status.divergence_of_velocity)
    du2.TW.current_time = tk1
    du2.TW.___DO_push_all_to_instant___()
    DIV_L2_error_tk1 = du2.error.L()
    KE2_tk1 = 0.5 * u2.DO.compute_L2_inner_product_energy_with(M=M2)
    E1_tk1  = 0.5 * w1.DO.compute_L2_inner_product_energy_with(M=M1)

    if tk == t0 and save_HS:
        P3.TW.current_time = t0 + dt / 2
        save(P3, f'pwc1_hsP3_N{N}k{k}elements{mesh.elements.GLOBAL_num}')



    LHS_HAL_00_A = LS_half_lhs00.assembled
    RHS_HAL_00_A = LS_half_B0.assembled
    LHS_HAL_00_A.DO.identify_rows(u1_WALL_DoFs)
    RHS_HAL_00_A.DO.clear_entries(u1_WALL_DoFs)
    LS_hal.lhs[0][0] = LHS_HAL_00_A
    LS_hal.rhs[0] = RHS_HAL_00_A
    X0 = LS_hal.solve.results
    LS_hal.solve('GMRES', 'MPI1')(X0, restart=restart, tol=tol, maxiter=maxiter, conserving=False)

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

    if abs(round(tk1) - tk1) < 1e-8: save([u1, u2, w1, w2], f'PWC1_u1u2w1w2_N{N}k{k}_t{int(round(tk1))}')

    u1.cochain.local = _u1_new_cochain_  # renew u1 cochain to time tk+half
    w2.cochain.local = u1.coboundary.cochain_local  # renew w2 cochain to time tk+half

    message = ['Inner solver: ' + LS_hal.solve.message, 'Outer solver: ' + LS_int.solve.message]

    return 1, 0, message, KE1_tk1, KE2_tk1, H1_tk1, H2_tk1, E1_tk1, E2_tk1, \
               u1u2_diff_tk1, w1w2_diff_tk1, DIV_L2_error_tk1



SI(SOLVER, [0, 0, 0, 0, 0, 0, 0, 0, None])
SI.run()

