# -*- coding: utf-8 -*-
"""
$ mpiexec python _3dCSCG\APP\contents\icpsNS\BC\test1.py

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from TOOLS.__DEPRECATED__.linear_system.main import LinearSystem
from TOOLS.linear_algebra.data_structures import GlobalMatrix, GlobalVector
from scipy import sparse as spspa

N = 1
k = 2
t0 = 0
t = 2
steps=100
nu=1
rho = 1


dt = (t-t0)/steps

quad_degree = [2 * N, 2 * N, 2 * N]
mesh = MeshGenerator('crazy', c=0.0, bounds=[(0, 1), (0, 1), (0, 1)])([k, k, k])
space = SpaceInvoker('polynomials')([('Lobatto', N), ('Lobatto', N), ('Lobatto', N)])
FC = FormCaller(mesh, space)
es = ExactSolutionSelector(mesh)('icpsNS:sincosCUCD1', nu=nu, rho=rho)

P0 = FC('0-f', is_hybrid=False, orientation='inner', name='inner-total-pressure',
        numbering_parameters={'scheme_name': 'Naive', })
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
mE01M1_A = - M1E10_A.T
mE01M1_A.DO.clear_row(-1)
lhs11_A = GlobalMatrix((P0.GLOBAL_num_dofs, P0.GLOBAL_num_dofs))
lhs11_A.DO.set_value(-1, -1, 1)


M1E10_A.DO.clear_row(0)
lhs00_Hf_A.DO.clear_row(0)
lhs00_Hf_A.DO.set_value(0, 0, 1)

lhs = [[lhs00_Hf_A, M1E10_A],
       [mE01M1_A  , lhs11_A]]

FI = f.___PRIVATE_do_inner_product_with_space_of___(u1, quad_degree=quad_degree)
FO = f.___PRIVATE_do_inner_product_with_space_of___(u2, quad_degree=quad_degree)

f.current_time = t0 + dt / 4
B0 = (2 * M1 / dt - 0.5 * CP1 - 0.5 * nu * E12M2E21) @ u1.cochain.EWC + FI
B0.gathering_matrix = u1
B0 = B0.assembled
B1 = GlobalVector(spspa.csc_matrix((P0.GLOBAL_num_dofs, 1)))
LS_hal = LinearSystem(lhs, rhs=[B0, B1])

# LS_hal.visualize.spy()

cn = LS_hal.condition.condition_number
if rAnk == mAster_rank:
    print('Inner system condition number:', cn)
    print(LS_hal.structure.regularity)