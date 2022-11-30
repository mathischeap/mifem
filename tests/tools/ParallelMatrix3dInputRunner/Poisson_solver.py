# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/19 9:26 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from __init__ import cscg2, tools
import numpy as np


def PoissonSolver(N, K, c, show_info=False):
       """

       Parameters
       ----------
       N :
       K :
       c :
       show_info :

       Returns
       -------
       p_error_L2 :
       u_error_L2 :

       """

       pB = ['Upper', 'Down', 'Left', 'Right']
       def p_fun(t, x, y): return  np.sin(np.pi * x) * np.sin(np.pi * y) + 0 * t
       def u_fun(t, x, y): return  np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) + 0 * t
       def v_fun(t, x, y): return  np.pi * np.sin(np.pi * x) * np.cos(np.pi * y) + 0 * t
       def F_fun(t, x, y): return 2 * np.pi ** 2 * np.sin(np.pi * x) * np.sin(np.pi * y) + 0 * t

       mesh = cscg2.mesh('crazy', c=c)([K, K], show_info=show_info)
       space = cscg2.space('polynomials')([('Lobatto', N), ('Lobatto', N)])
       fc = cscg2.form(mesh, space)

       p = fc('2-f-o', name='potential', is_hybrid=True)
       u = fc('1-f-o', name='velocity', is_hybrid=True)
       f = fc('2-f-o', name='source', is_hybrid=True)
       t = fc('1-t-o', name='trace')

       t.BC.boundaries = pB

       s = fc('scalar', F_fun)
       f.CF = s
       f.CF.current_time = 0
       f.discretize()

       m2 = p.matrices.mass
       m1 = u.matrices.mass
       e21 = u.matrices.incidence
       e12 = e21.T
       n = t.matrices.trace

       e12m2 = e12 @ m2
       a = tools.ewc.bmat((
              [m1, e12m2, -n.T],
              [e21, None, None],
              [n, None, None]))
       a.assembler.chain_method = 'sequent'
       a.gathering_matrices = ((u, p, t), (u, p, t))
       b0 = tools.ewc.vector(u)
       b1 = -f.cochain.EWC
       b2 = tools.ewc.vector(t)
       b = tools.ewc.concatenate([b0, b1, b2])
       b.assembler.chain_method = 'sequent'
       b.gathering_matrix = (u, p, t)
       ls = tools.milinalg.LinearSystem(a, b)
       a.customize.identify_global_rows_according_to(2, t.BC.interpret)
       results = ls.solve('direct')()[0]
       results.do.distributed_to(u, p, t, chain_method=b.gathering_matrix)

       ps = fc('scalar', p_fun)
       uv = fc('vector', [u_fun, v_fun])
       p.CF = ps
       p.CF.current_time = 0
       p_error_L2 = p.error.L()
       u.CF = uv
       u.CF.current_time = 0
       u_error_L2 = u.error.L()
       return p_error_L2, u_error_L2


if __name__ == "__main__":
    # mpiexec -n 4 python tests/tools/ParallelMatrix3dInputRunner/Poisson_solver.py
    R = PoissonSolver(3, 5, 0)
    print(R)