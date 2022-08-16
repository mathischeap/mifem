# -*- coding: utf-8 -*-
"""Use the Schur solver to solve the hybrid problem.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/15 8:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from tools.deprecated.linear_system.main import LinearSystem
from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector, GlobalMatrix
from scipy import sparse as spspa


def FUNC_scalar_Laplace_solver(N, K, c, show_info=True):
       """

       Parameters
       ----------
       N :
       K :
       c :
       show_info : bool

       Returns
       -------
       u_error_L2 :
       p_error_L2 :
       du_error_L2 :

       """
       mesh = MeshGenerator('crazy', c=c, bounds=[(0,1), (0,1)])([K, K], EDM='debug', show_info=show_info)
       space = SpaceInvoker('polynomials')([N, N], show_info=show_info)
       FC = FormCaller(mesh, space)
       ES = ExactSolutionSelector(mesh)('sL:sincos1')

       u = FC('1-f-o', is_hybrid=False)
       p = FC('2-f-o', is_hybrid=False)
       f = FC('2-f-o', is_hybrid=False)

       B0 = GlobalVector(spspa.csc_matrix((u.num.GLOBAL_dofs, 1)))
       f.TW.func.___DO_set_func_body_as___(ES, "source_term")
       f.TW.___DO_push_all_to_instant___(0)
       f.discretize()
       B1 = f.cochain.EWC
       B1.gathering_matrix = p
       B1 = B1.assembled

       M1 = u.matrices.mass
       M2 = p.matrices.mass
       E21 = u.matrices.incidence
       E12 = E21.T
       E12M2 = E12 @ M2

       M1.gathering_matrices = (u, u)
       E21.gathering_matrices = (p, u)
       E12M2.gathering_matrices = (u, p)

       M1 = M1.assembled
       E21 = E21.assembled
       E12M2 = E12M2.assembled
       lhs22 = GlobalMatrix((p.num.GLOBAL_dofs, p.num.GLOBAL_dofs))

       lhs = ([M1, E12M2 ],
              [E21, lhs22])

       LS = LinearSystem(lhs, rhs=[B0, B1])

       LS.solve('Direct', 'spspalinalg')()

       LS.solve.results.do.distributed_to(u, p)

       u.TW.func.___DO_set_func_body_as___(ES, "velocity")
       u.TW.___DO_push_all_to_instant___(0)
       u_error_L2 = u.error.L()

       p.TW.func.___DO_set_func_body_as___(ES, "potential")
       p.TW.___DO_push_all_to_instant___(0)
       p_error_L2 = p.error.L()

       du_error_L2 = u.do.compute_Ln_norm_of_coboundary()

       return u_error_L2, p_error_L2, du_error_L2

if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/auxiliaries/SL_eBC.py
    R = FUNC_scalar_Laplace_solver(3, 5, 0)
    print(R)