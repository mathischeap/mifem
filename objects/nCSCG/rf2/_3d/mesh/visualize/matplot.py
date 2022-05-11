# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 12:40 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.freeze.base import FrozenOnly
import matplotlib.pyplot as plt
import matplotlib


class _3nCSCG_MeshVisualizeMatplot(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, saveto=None, aspect='equal',
                 usetex=False,
                 labelsize=15, ticksize=15):
        """"""
        CPD = dict()
        for ind in self._mesh_:
            cell = self._mesh_(ind)
            assert cell.___isroot___
            CPD[str(ind)] = cell.coordinate_transformation.___PRIVATE_plot_data___()

        CPD = cOmm.gather(CPD, root=mAster_rank)

        if rAnk != mAster_rank: return

        ___ = dict()
        for cpd in CPD:
            ___.update(cpd)
        CPD = ___

        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for ind in CPD:
            lines = CPD[ind]
            for xyz in lines:
                plt.plot(*xyz, color='k', linewidth=0.75)

                if aspect == 'equal':
                    x, y, z = xyz
                    x_lim.append(np.max(x))
                    x_lim.append(np.min(x))
                    y_lim.append(np.max(y))
                    y_lim.append(np.min(y))
                    z_lim.append(np.max(z))
                    z_lim.append(np.min(z))

        if aspect == 'equal': ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))
        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        #---------------------- save to ---------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')
        plt.close()







if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rf2/_3d/mesh/visualize/matplot.py
    from objects.nCSCG.rf2._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None)
    if 0 in mesh.cscg.elements:
        c0 = mesh(0)
        c0.do.refine()

    mesh.visualize()