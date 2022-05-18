# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 12:40 PM
"""
import sys

import numpy as np

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

    def __call__(self, density=None, color_space=True, show_indices=False,
        saveto=None, usetex=False, labelsize=12, ticksize=12):
        """

        Parameters
        ----------
        density
        color_space
        show_indices
        saveto
        usetex
        labelsize
        ticksize

        Returns
        -------

        """
        CPD = dict()
        indices_dict = dict()
        color_Space_dict = dict()
        for ci in self._mesh_: # ao through all local cell indices.
            cell = self._mesh_(ci)
            assert cell.___isroot___
            CPD[repr(cell)] = cell.coordinate_transformation.___PRIVATE_plot_data___(density=density)

            if show_indices:
                indices = cell.indices
                center_coo = cell.coordinate_transformation.mapping(0, 0)
                indices_dict[repr(cell)] = [indices, center_coo]

            if color_space:
                color_Space_dict[repr(cell)] = cell.space.N

        CPD = cOmm.gather(CPD, root=mAster_rank)
        if show_indices:
            indices_dict = cOmm.gather(indices_dict, root=mAster_rank)
        if color_space:
            color_Space_dict = cOmm.gather(color_Space_dict, root=mAster_rank)


        if rAnk != mAster_rank: return

        ___ = dict()
        for _ in CPD:
            ___.update(_)
        CPD = ___

        if show_indices:
            ___ = dict()
            for _ in indices_dict:
                ___.update(_)
            indices_dict = ___

        if color_space:
            ___ = dict()
            for _ in color_Space_dict:
                ___.update(_)
            color_Space_dict = ___

            Ns = set(color_Space_dict.values())
            Max_N = max(Ns)
            Min_N = min(Ns)
            if Max_N == Min_N:
                homogeneous_N = True
            else:
                homogeneous_N = False
                dN = Max_N - Min_N

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
            for xy in lines:
                plt.plot(*xy, color='k', linewidth=0.5)

            if show_indices:
                text, coo = indices_dict[ind]
                plt.text(*coo, text, ha='center', va='center')

            if color_space:
                if homogeneous_N:
                    pass
                else:
                    X, Y = [None, None, None, None], [None, None, None, None]
                    for xy, i in zip(lines, [0, 2, 3, 1]):
                        x, y = xy
                        if i == 0:
                            X[i] = x
                            Y[i] = y
                        elif i == 1:
                            X[i] = x[1:]
                            Y[i] = y[1:]
                        elif i == 2:
                            X[i] = x[::-1][1:]
                            Y[i] = y[::-1][1:]
                        else:
                            X[i] = x[::-1][1:]
                            Y[i] = y[::-1][1:]
                    X = np.concatenate(X)
                    Y = np.concatenate(Y)

                    N = color_Space_dict[ind]
                    c = 0.9 - (N - Min_N) * 0.8 / dN
                    plt.fill(X, Y, color=(c, c, c, 0.3))

        #---------------------- save to --------------------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')
        plt.close()






if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rfT2/_2d/mesh/visualize/matplot.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100)
    # from objects.nCSCG.rfT2._2d.master import MeshGenerator
    # mesh = MeshGenerator('crazy', c=0.)([3, 3], 2, EDM='chaotic', show_info=False)
    # mesh.do.unlock()
    # i = 4
    # if i in mesh.cscg.elements:
    #     c = mesh(i)
    #     c.do.refine()
    #     c0 = c(0)
    #     c0.do.refine()
    #     c3 = c(3)
    #     c3.do.refine()

    mesh.visualize(show_indices=False)