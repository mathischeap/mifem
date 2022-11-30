# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:25 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from __init__ import miTri, tools
import numpy as np

from components.freeze.base import FrozenOnly
from components.miscellaneous.miprint import miprint
from components.miscellaneous.mios import remove, isfile
from tools.run.reader import ParallelMatrix3dInputRunner

import os
current_dir = os.path.dirname(__file__)

def __PoissonSolver1__(N, K, c):
       """

       Parameters
       ----------
       N :
       K :
       c :

       Returns
       -------
       p_error_L2 :
       u_error_L2 :
       u_error_Hdiv :
       mass_conservation :

       """

       def p_fun(t, x, y): return np.sin(np.pi * x) * np.sin(np.pi * y) + 0 * t * c
       def u_fun(t, x, y): return np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) + 0 * t
       def v_fun(t, x, y): return np.pi * np.sin(np.pi * x) * np.cos(np.pi * y) + 0 * t
       def F_fun(t, x, y): return 2 * np.pi ** 2 * np.sin(np.pi * x) * np.sin(np.pi * y) + 0 * t

       fc = miTri.call(f'st{K}', N)

       p = fc('2-f-o', name='potential')
       u = fc('1-f-o', name='velocity')
       f = fc('2-f-o', name='source')

       s = fc('scalar', F_fun)
       s.current_time = 0
       f.CF = s
       f.discretize()

       m2 = p.matrices.mass
       m1 = u.matrices.mass

       e21 = u.matrices.incidence
       e12 = e21.T

       e12m2 = e12 @ m2
       a = tools.ewc.bmat(([m1 , e12m2],
                              [e21, None ]))

       a.assembler.chain_method = 'sequent'
       a.gathering_matrices = ((u, p), (u, p))

       b0 = tools.ewc.vector(fc.mesh.elements, u.num.basis)
       b1 = -f.cochain.EWC
       b = tools.ewc.concatenate([b0, b1])
       b.assembler.chain_method = 'sequent'
       b.gathering_matrix = (u, p)
       ls = tools.milinalg.LinearSystem(a, b)

       results = ls.solve('direct')()[0]
       results.do.distributed_to(u, p, chain_method='sequent')

       ps = fc('scalar', p_fun)
       uv = fc('vector', [u_fun, v_fun])
       ps.current_time = 0
       uv.current_time = 0

       p.CF = ps
       p_error_L2 = p.error.L()
       u.CF = uv
       u_error_L2 = u.error.L()
       u_error_Hdiv = u.error.H()

       du = u.coboundary()
       du_plus_f = du + f
       mass_conservation = du_plus_f.do.compute_Ln_norm()
       assert mass_conservation < 1e-11
       if K >= 10 and N == 3:
           assert p_error_L2 < 1e-3
           assert u_error_L2 < 1e-3

       return p_error_L2, u_error_L2, u_error_Hdiv, mass_conservation



class miUsGrid_Triangle_Poisson(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("miUsTri [miUsGrid_Triangle_Poisson] ....", flush=True)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        Ns = [[1, 1, 1, 1, 1,],
              [2, 2, 2, 2, 2,],
              [3, 3, 3, 3, 3,]]
        Ks = [[2, 4, 6, 8, 10,],
              [2, 4, 6, 8, 10,],
              [2, 4, 6, 8, 10,], ]
        cs = [0, ]

        pr = ParallelMatrix3dInputRunner(__PoissonSolver1__)

        if isfile(current_dir + '/TriangularPoissonTest.txt'): remove(current_dir + '/TriangularPoissonTest.txt')

        pr.iterate(Ns, Ks, cs, writeto=current_dir + '/TriangularPoissonTest.txt', show_progress=False)

        remove(current_dir + '/TriangularPoissonTest.txt')

        pr.visualize('loglog', 'N', 'p_error_L2', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \varphi^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 0},
                                          1: {'tp': (0.02, 0.2), 'order': 1},
                                          2: {'tp': (0.02, 0.2), 'order': 2}},
                     saveto=current_dir + '/p_error_L2.png')

        pr.visualize('loglog', 'N', 'u_error_L2', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \boldsymbol{u}^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                          1: {'tp': (0.02, 0.2), 'order': 2},
                                          2: {'tp': (0.02, 0.2), 'order': 3}},
                     saveto=current_dir + '/u_error_L2.png')

        pr.visualize('loglog', 'N', 'u_error_Hdiv', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \boldsymbol{u}^h\right\|_{H(\mathrm{div})-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 0},
                                          1: {'tp': (0.02, 0.2), 'order': 1},
                                          2: {'tp': (0.02, 0.2), 'order': 2}},
                     saveto=current_dir + '/u_error_Hdiv.png')

        pr.visualize('loglog', 'N', 'mass_conservation', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \mathrm{d}\boldsymbol{u}^h + f^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/mass_conservation_L2.png')

        a, b = pr.results.to_numpy()[:,3][-2:]
        a = np.log10(a)
        b = np.log10(b)
        c = np.log10(1/8)
        d = np.log10(1/10)

        assert abs((a-b)/(c-d) - 2) < 0.05, f"[N-1]th order convergence rate lost."

        return 1



if __name__ == "__main__":
    # mpiexec -n 4 python tests/objects/miUsGrid/triangular/unittests/standardForms/Poisson/test.py
    miUsGrid_Triangle_Poisson()()
