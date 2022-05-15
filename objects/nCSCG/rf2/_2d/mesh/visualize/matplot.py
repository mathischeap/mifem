# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 12:40 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm
from screws.freeze.base import FrozenOnly
import matplotlib.pyplot as plt
import matplotlib




class _2nCSCG_MeshVisualizeMatplot(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, show_indices=False,
        saveto=None, usetex=False, labelsize=12, ticksize=12):
        """"""
        CPD = dict()
        indices_dict = dict()
        for ci in self._mesh_: # ao through all local cell indices.
            cell = self._mesh_(ci)
            assert cell.___isroot___
            CPD[str(ci)] = cell.coordinate_transformation.___PRIVATE_plot_data___()

            if show_indices:
                indices = cell.indices
                center_coo = cell.coordinate_transformation.mapping(0, 0)
                indices_dict[str(ci)] = [indices, center_coo]

        CPD = cOmm.gather(CPD, root=mAster_rank)
        indices_dict = cOmm.gather(indices_dict, root=mAster_rank)

        if rAnk != mAster_rank: return

        ___ = dict()
        for _ in CPD:
            ___.update(_)
        CPD = ___

        ___ = dict()
        for _ in indices_dict:
            ___.update(_)
        indices_dict = ___

        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        plt.xlabel(r"$x$", fontsize=labelsize)
        plt.ylabel(r"$y$", fontsize=labelsize)
        plt.tick_params(axis='both', which='both', labelsize=ticksize)

        for ind in CPD:
            lines = CPD[ind]
            for xyz in lines:
                plt.plot(*xyz, color='k', linewidth=0.5)
            if show_indices:
                text, coo = indices_dict[ind]
                plt.text(*coo, text, ha='center', va='center')

        #---------------------- save to --------------------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')
        plt.close()






if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rf2/_2d/mesh/visualize/matplot.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy', c=0.)([5, 5], 2, EDM='chaotic', show_info=False)
    mesh.do.unlock()

    i = 4
    if i in mesh.cscg.elements:
        c = mesh(i)
        c.do.refine()
        c0 = c(0)
        c0.do.refine()
        c3 = c(3)
        c3.do.refine()

    mesh.visualize()