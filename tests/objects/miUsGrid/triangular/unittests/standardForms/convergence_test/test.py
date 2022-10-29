# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/26 7:47 PM
"""
import os
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly
import numpy as np

from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mios import remove, isfile
from tools.run.reader import ParallelMatrix3dInputRunner
from objects.miUsGrid.triangular.master import Call

current_dir = os.path.dirname(__file__)

def func(t, x, y): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + t
def fx(t, x, y): return np.sin(2*np.pi*x) * np.cos(np.pi*y) + t
def fy(t, x, y): return np.cos(2*np.pi*x) * np.sin(np.pi*y) + t


# noinspection PyUnusedLocal
def ___Pr_error_function___(N, K, c):
    """

    Parameters
    ----------
    N :
    K :
    c :

    Returns
    -------
    f0i_er :
    f0o_er :
    f1i_er :
    f1o_er :
    f2i_er :
    f2o_er :
    df0i_er :
    df0o_er :
    df1i_er :
    df1o_er :

    """
    mesh_id = f"st{int(K)}"
    fc = Call(mesh_id, N)
    if K == 8:
        fc.mesh.visualize(saveto=current_dir + '/hc_mesh_K8.pdf', show_singular_vertex=False)

    f0i = fc('0-f-i')
    f0o = fc('0-f-o')
    f1i = fc('1-f-i')
    f1o = fc('1-f-o')
    f2i = fc('2-f-i')
    f2o = fc('2-f-o')

    scalar = fc('scalar', func)
    vector = fc('vector', (fx, fy))

    scalar.current_time = 0
    vector.current_time = 0

    f0i.CF = scalar
    f0o.CF = scalar
    f1i.CF = vector
    f1o.CF = vector
    f2i.CF = scalar
    f2o.CF = scalar

    f0i.discretize()
    f0o.discretize()
    f1i.discretize()
    f1o.discretize()
    f2i.discretize()
    f2o.discretize()

    f0i_er = f0i.error.L()
    f0o_er = f0o.error.L()
    f1i_er = f1i.error.L()
    f1o_er = f1o.error.L()
    f2i_er = f2i.error.L()
    f2o_er = f2o.error.L()

    grad_scalar = scalar.numerical.grad
    curl_scalar = scalar.numerical.curl

    div_vector = vector.numerical.div
    rot_vector = vector.numerical.rot

    grad_scalar.current_time = 0
    curl_scalar.current_time = 0
    div_vector.current_time = 0
    rot_vector.current_time = 0

    df0i = f0i.coboundary()
    df0i.CF = grad_scalar
    df0i_er = df0i.error.L()

    df0o = f0o.coboundary()
    df0o.CF = curl_scalar
    df0o_er = df0o.error.L()

    df1i = f1i.coboundary()
    df1i.CF = rot_vector
    df1i_er = df1i.error.L()

    df1o = f1o.coboundary()
    df1o.CF = div_vector
    df1o_er = df1o.error.L()

    return f0i_er, f0o_er, f1i_er, f1o_er, f2i_er, f2o_er, df0i_er, df0o_er, df1i_er, df1o_er

# noinspection PyUnusedLocal
def ___Pr_error_function_pc___(N, K, c):
    """

    Parameters
    ----------
    N :
    K :
    c :

    Returns
    -------
    f0i_er :
    f0o_er :
    f1i_er :
    f1o_er :
    f2i_er :
    f2o_er :
    df0i_er :
    df0o_er :
    df1i_er :
    df1o_er :

    """
    fc = Call('rand0', N)

    if N == 10: fc.mesh.visualize(saveto=current_dir + '/pc_mesh.pdf')

    f0i = fc('0-f-i')
    f0o = fc('0-f-o')
    f1i = fc('1-f-i')
    f1o = fc('1-f-o')
    f2i = fc('2-f-i')
    f2o = fc('2-f-o')

    scalar = fc('scalar', func)
    vector = fc('vector', (fx, fy))

    scalar.current_time = 0
    vector.current_time = 0

    f0i.CF = scalar
    f0o.CF = scalar
    f1i.CF = vector
    f1o.CF = vector
    f2i.CF = scalar
    f2o.CF = scalar

    f0i.discretize()
    f0o.discretize()
    f1i.discretize()
    f1o.discretize()
    f2i.discretize()
    f2o.discretize()

    f0i_er = f0i.error.L()
    f0o_er = f0o.error.L()
    f1i_er = f1i.error.L()
    f1o_er = f1o.error.L()
    f2i_er = f2i.error.L()
    f2o_er = f2o.error.L()

    grad_scalar = scalar.numerical.grad
    curl_scalar = scalar.numerical.curl

    div_vector = vector.numerical.div
    rot_vector = vector.numerical.rot

    grad_scalar.current_time = 0
    curl_scalar.current_time = 0
    div_vector.current_time = 0
    rot_vector.current_time = 0

    df0i = f0i.coboundary()
    df0i.CF = grad_scalar
    df0i_er = df0i.error.L()

    df0o = f0o.coboundary()
    df0o.CF = curl_scalar
    df0o_er = df0o.error.L()

    df1i = f1i.coboundary()
    df1i.CF = rot_vector
    df1i_er = df1i.error.L()

    df1o = f1o.coboundary()
    df1o.CF = div_vector
    df1o_er = df1o.error.L()

    if N >= 20: # a check of reduction and reconstruction.
        # when N is significantly large, we should reach machine zero.
        assert f0i_er < 1e-11
        assert f0o_er < 1e-11
        assert f1i_er < 1e-11
        assert f1o_er < 1e-11
        assert f2i_er < 1e-11
        assert f2o_er < 1e-11
        assert df1i_er < 1e-8
        assert df1o_er < 1e-8
        assert df0i_er < 1e-8
        assert df0o_er < 1e-8

    return f0i_er, f0o_er, f1i_er, f1o_er, f2i_er, f2o_er, df0i_er, df0o_er, df1i_er, df1o_er


class miUsGrid_TriangleMesh_ConvergenceTest(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("miUsTri [miUsGrid_TriangleMesh_ConvergenceTest] ....", flush=True)
        self._freeze_self_()


    def __call__(self):
        """"""

        Ns = [[1, 1, 1, 1, 1,],
              [2, 2, 2, 2, 2,],
              [3, 3, 3, 3, 3,]]
        Ks = [[2, 4, 6, 8, 10,],
              [2, 4, 6, 8, 10,],
              [2, 4, 6, 8, 10,],]
        cs = [0,]

        pr = ParallelMatrix3dInputRunner(___Pr_error_function___)

        if isfile(current_dir + '/WTP.txt'): remove(current_dir + '/WTP.txt')

        pr.iterate(Ns, Ks, cs, writeto=current_dir + '/WTP.txt', show_progress=False)



        remove(current_dir + '/WTP.txt')

        pr.visualize('loglog', 'N', 'f0i_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \varphi^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 2},
                                          1: {'tp': (0.02, 0.2), 'order': 3},
                                          2: {'tp': (0.02, 0.2), 'order': 4}},
                     saveto=current_dir + '/f0_error_L2.pdf')

        pr.visualize('loglog', 'N', 'f1i_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \boldsymbol{\omega}^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                          1: {'tp': (0.02, 0.2), 'order': 2},
                                          2: {'tp': (0.02, 0.2), 'order': 3}},
                     saveto=current_dir + '/f1i_error_L2.pdf')

        pr.visualize('loglog', 'N', 'f1o_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \boldsymbol{u}^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                          1: {'tp': (0.02, 0.2), 'order': 2},
                                          2: {'tp': (0.02, 0.2), 'order': 3}},
                     saveto=current_dir + '/f1o_error_L2.pdf')

        pr.visualize('loglog', 'N', 'f2i_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| f^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 0},
                                          1: {'tp': (0.02, 0.2), 'order': 1},
                                          2: {'tp': (0.02, 0.2), 'order': 2}},
                     saveto=current_dir + '/f2_error_L2.pdf')

        pr.visualize('loglog', 'N', 'df0i_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \mathrm{grad}(\varphi^h)\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                          1: {'tp': (0.02, 0.2), 'order': 2},
                                          2: {'tp': (0.02, 0.2), 'order': 3}},
                     saveto=current_dir + '/df0i_error_L2.pdf')

        pr.visualize('loglog', 'N', 'df0o_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \mathrm{curl}(\varphi^h)\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                          1: {'tp': (0.02, 0.2), 'order': 2},
                                          2: {'tp': (0.02, 0.2), 'order': 3}},
                     saveto=current_dir + '/df0o_error_L2.pdf')

        pr.visualize('loglog', 'N', 'df1i_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \mathrm{rot}(\boldsymbol{\omega})^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 0},
                                          1: {'tp': (0.02, 0.2), 'order': 1},
                                          2: {'tp': (0.02, 0.2), 'order': 2}},
                     saveto=current_dir + '/df1i_error_L2.pdf')

        pr.visualize('loglog', 'N', 'df1o_er', prime='input2', hcp=1, usetex=True,
                     labels=['$N=1$', '$N=2$', '$N=3$',],
                     styles=["-s", "-v", '-^'],
                     colors=[(0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                     title=False,
                     # yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3],
                     xlabel=r'$1/K$',
                     ylabel=r"$\left\| \mathrm{div}(\boldsymbol{u})^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 0},
                                          1: {'tp': (0.02, 0.2), 'order': 1},
                                          2: {'tp': (0.02, 0.2), 'order': 2}},
                     saveto=current_dir + '/df1o_error_L2.pdf')

        Ns = [[i for i in range(1, 18)],]
        Ks = [[0 for _ in range(1, 18)],]
        cs = [0, ]
        pr = ParallelMatrix3dInputRunner(___Pr_error_function_pc___)

        if isfile(current_dir + '/WTP_pc.txt'): remove(current_dir + '/WTP_pc.txt')

        pr.iterate(Ns, Ks, cs, writeto=current_dir + '/WTP_pc.txt', show_progress=False)
        remove(current_dir + '/WTP_pc.txt')

        pr.visualize('semilogy', 'K', 'f0i_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \varphi^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_f0_error_L2.pdf')

        pr.visualize('semilogy', 'K', 'f1i_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \boldsymbol{\omega}^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_f1i_error_L2.pdf')

        pr.visualize('semilogy', 'K', 'f1o_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \boldsymbol{u}^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_f1o_error_L2.pdf')

        pr.visualize('semilogy', 'K', 'f2i_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| f^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_f2_error_L2.pdf')


        pr.visualize('semilogy', 'K', 'df0i_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \mathrm{grad}(\varphi)^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_df0i_error_L2.pdf')

        pr.visualize('semilogy', 'K', 'df0o_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \mathrm{curl}(\varphi)^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_df0o_error_L2.pdf')


        pr.visualize('semilogy', 'K', 'df1i_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \mathrm{rot}(\boldsymbol{\omega})^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_df1i_error_L2.pdf')

        pr.visualize('semilogy', 'K', 'df1o_er', prime='input2', usetex=True,
                     labels=False,
                     styles=["-v", "-v"],
                     colors=[(0.4, 0.4, 0.4, 1),], title=False,
                     xlabel=r'$N$',
                     ylabel=r"$\left\| \mathrm{div}(\boldsymbol{u})^h\right\|_{L^2-\mathrm{error}}$",
                     order_text_size=15,
                     saveto=current_dir + '/pc_df1o_error_L2.pdf')

        return 1






if __name__ == "__main__":
    # mpiexec -n 4 python tests/objects/miUsGrid/triangular/unittests/standard_forms/convergence_test/test.py
    miUsGrid_TriangleMesh_ConvergenceTest()()
