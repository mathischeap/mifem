# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/15 8:17 PM
"""
import os
import sys

if './' not in sys.path: sys.path.append('./')
from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mios import remove

from objects.CSCG._2d.__tests__.unittests.auxiliaries.SL_eBC import FUNC_scalar_Laplace_solver
from tools.run.reader import ParallelMatrix3dInputRunner, RunnerDataReader

def WellTest_ParallelMatrix3dInputRunner():
    miprint("WTP [WellTest_ParallelMatrix3dInputRunner] ....", flush=True)
    current_dir = os.path.dirname(__file__)

    Ns = [[1, 1, 1, 1, 1],
          [3, 3, 3, 3, 3]]
    Ks = [[2, 4, 6, 8, 10],
          [2, 4, 6, 8, 10]]
    cs = [0, 0.1, 0.2]

    pr = ParallelMatrix3dInputRunner(FUNC_scalar_Laplace_solver)
    # pr.iterate(Ns, Ks, cs, writeto='WTP.txt', show_info=False, show_progress=False)
    pr.iterate(Ns, Ks, cs, writeto=current_dir + '/WTP.txt', show_info=False, show_progress=False)
    # PR = RunnerDataReader('WTP.txt')
    PR = RunnerDataReader(current_dir + '/WTP.txt')
    PR.visualize.quick('N', y='u_error_L2', saveto=current_dir + '/u_error_L2_quick.png')
    PR.visualize.quick('K', y='p_error_L2', saveto=current_dir + '/p_error_L2_quick.png')
    #
    PR.visualize('loglog', 'N', 'u_error_L2', prime='input2', hcp=1, usetex=True,
             labels=['$N=1,c=0$'  , '$N=3, c=0$',
                    '$N=1,c=0.1$', '$N=3, c=0.1$',
                    '$N=1,c=0.2$', '$N=3, c=0.2$'],
             styles=["-s", "-v"],
             colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
             title=False,
             yticks=[1e1, 1e0,  1e-1, 1e-2, 1e-3],
             xlabel=r'$1/K$',
             ylabel=r"$\left\| \boldsymbol{u}_h\right\|_{L^2-\mathrm{error}}$",
             order_text_size=15,
             plot_order_triangle={0: {'tp':(0.02,-0.5),'order':1},
                                 1: {'tp':(0.02, 0.2),'order':3}},
             saveto=current_dir + '/u_error_L2.png')
    PR.visualize('loglog', 'N', 'p_error_L2', prime='input2', hcp=1, usetex=True,
                 labels=['$N=1,c=0$' , '$N=3, c=0$',
                         '$N=1,c=0.1$', '$N=3, c=0.1$',
                         '$N=1,c=0.2$', '$N=3, c=0.2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 yticks=[1e0, 1e-1, 1e-2, 1e-3, 1e-4],
                 xlabel=r'$1/K$',
                 ylabel=r"$\left\| p_h\right\|_{L^2-\mathrm{error}}$", order_text_size=15,
                 plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                      1: {'tp': (0.02, 0.2), 'order': 3}},
                 saveto=current_dir + '/p_error_L2.png')

    PR.visualize('loglog', 'N', 'du_error_L2', prime='input2', hcp=1, usetex=True,
                 labels=['$N=1,c=0$' , '$N=3, c=0$',
                         '$N=1,c=0.1$', '$N=3, c=0.1$',
                         '$N=1,c=0.2$', '$N=3, c=0.2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 xlabel=r'$1/K$',
                 ylabel=r"$\left\| \mathrm{d}\boldsymbol{u}_h\right\|_{L^2-\mathrm{error}}$",
                 saveto=current_dir + '/du_error_L2.png')



    Ns = [[1, 2, 3, 4, 5],
          [1, 2, 3, 4, 5]]
    Ks = [[2, 2, 2, 2, 2],
          [4, 4, 4, 4, 4]]
    cs = [0, 0.1, 0.2]
    pr = ParallelMatrix3dInputRunner(FUNC_scalar_Laplace_solver)
    pr.iterate(Ns, Ks, cs, writeto=current_dir + '/WTP_Pc.txt', show_info=False, show_progress=False)
    pr.visualize('semilogy', 'K', 'u_error_L2', prime='input2', usetex=True,
                 labels=['$K=2,c=0$'  , '$K=4, c=0$',
                         '$K=2,c=0.1$', '$K=4, c=0.1$',
                         '$K=2,c=0.2$', '$K=4, c=0.2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)], title=False,
                 yticks=[1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4],
                 xlabel=r'$N$',
                 ylabel=r"$\left\| \boldsymbol{u}_h\right\|_{L^2-\mathrm{error}}$", order_text_size=15,
                 saveto=current_dir + '/u_error_L2_Pc.png')

    pr.visualize('semilogy', 'K', 'p_error_L2', prime='input2', usetex=True,
                 labels=['$K=2,c=0$'  , '$K=4, c=0$',
                         '$K=2,c=0.1$', '$K=4, c=0.1$',
                         '$K=2,c=0.2$', '$K=4, c=0.2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)], title=False,
                 yticks=[1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5],
                 xlabel=r'$N$',
                 ylabel=r"$\left\| \boldsymbol{u}_h\right\|_{L^2-\mathrm{error}}$", order_text_size=15,
                 saveto=current_dir + '/p_error_L2_Pc.png')

    pr.visualize('semilogy', 'K', 'du_error_L2', prime='input2', usetex=True,
                 labels=['$K=2,c=0$'  , '$K=4, c=0$',
                         '$K=2,c=0.1$', '$K=4, c=0.1$',
                         '$K=2,c=0.2$', '$K=4, c=0.2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)], title=False,
                 xlabel=r'$N$',
                 ylabel=r"$\left\| \mathrm{d}\boldsymbol{u}_h\right\|_{L^2-\mathrm{error}}$", order_text_size=15,
                 saveto=current_dir + '/du_error_L2_Pc.png')

    remove(current_dir + '/WTP.txt')
    remove(current_dir + '/WTP_Pc.txt')

    return 1











if __name__ == "__main__":
    # mpiexec -n 4 python __tests__/unittests/Parallel_Matrix3dInput_Runner/WellTest.py
    WellTest_ParallelMatrix3dInputRunner()
