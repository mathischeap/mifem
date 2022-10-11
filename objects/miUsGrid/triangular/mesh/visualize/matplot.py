# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 3:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, mAster_rank, cOmm
import matplotlib
import matplotlib.pyplot as plt


class miUsGrid_TriangularMesh_Matplot(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.grid(*args, **kwargs)

    def grid(self, usetex=False, saveto=None, dpi=210, axis_on=True, show_singular_vertex=True):
        """"""
        grid_data = self._mesh_.elements.do.generate_grid_data()
        grid_data = cOmm.gather(grid_data, root=mAster_rank)
        if rAnk != mAster_rank: return
        GD = dict()
        for gd in grid_data:
            GD.update(gd)
        grid_data = GD

        #------- config plt ----------------------------------------------------------------------1
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        if axis_on:
            plt.tick_params(which='both', labeltop=False, labelright=False, top=False, right=False)
            plt.tick_params(axis='both', which='minor', direction='out', length=4)
            plt.tick_params(axis='both', which='major', direction='out', length=8)
            plt.tick_params(axis='both', which='both', labelsize=12)
            plt.tick_params(axis='x', which='both', pad=4)
            plt.tick_params(axis='y', which='both', pad=4)
            plt.xlabel('$x$', fontsize=12)
            plt.ylabel('$y$', fontsize=12)
        else:
            ax.set_axis_off()

        #------------------ grid data ------------------------------------------------------------1
        for i in grid_data:
            edges_coordinates, center, singular_vertex, inner_edges = grid_data[i]
            plt.plot(*edges_coordinates, '-k', linewidth=1)
            if show_singular_vertex:
                cx, cy = center
                sx, sy = singular_vertex
                plt.plot([cx, sx], [cy, sy], '--r', linewidth=0.5)
                plt.text(*center, f"{i}", color='gray', fontsize=11, ha='center', va='center')

                U, D, L, R = inner_edges
                plt.plot(*D, '--r', linewidth=0.5)
                plt.plot(*U, '--g', linewidth=0.5)
                plt.plot(*R, '--b', linewidth=0.5)
                plt.plot(*L, '--k', linewidth=0.5)

        #---------------------- save to ----------------------------------------------------------1
        if saveto is None:
            plt.show()
        else:
            if saveto[-4:] == '.pdf':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')
        plt.close()



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
