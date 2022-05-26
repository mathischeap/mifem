# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 8:28 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.freeze.base import FrozenOnly
import matplotlib.pyplot as plt
import matplotlib


class mpRfT2_Mesh_FutureRefinements_Preview(FrozenOnly):
    """"""

    def __init__(self, future):
        """"""
        self._future_ = future
        self._freeze_self_()

    def __call__(self, density=None, color_space=True, show_indices=False, show_boundaries=True,
        saveto=None, usetex=False, labelsize=12, ticksize=12):
        """

        Parameters
        ----------
        density
        color_space
        show_indices
        show_boundaries
        saveto
        usetex
        labelsize
        ticksize

        Returns
        -------

        """
        mesh = self._future_._mesh_
        cscg = mesh.cscg
        dN = mesh.dN
        ghost_cells = self._future_._mesh_.___Pr_initialize_cells___(dN, self._future_.rfd)

        CPD = dict()
        indices_dict = dict()
        color_Space_dict = dict()

        for i in ghost_cells:
            cell = ghost_cells[i]
            for j in cell:
                if len(j) == 1: # root-ghost-cell
                    ghost_root_cell = ghost_cells[j[0]]
                else:
                    ghost_basic_cell = ghost_cells[j[0]]
                    ghost_root_cell = ghost_basic_cell[j[1:]]

                assert ghost_root_cell.___isroot___
                CPD[repr(ghost_root_cell)] = ghost_root_cell.\
                    coordinate_transformation.___PRIVATE_plot_data___(density=density)

                if show_indices:
                    indices = ghost_root_cell.indices
                    center_coo = ghost_root_cell.coordinate_transformation.mapping(0, 0)
                    indices_dict[repr(ghost_root_cell)] = [indices, center_coo]

                if color_space:
                    color_Space_dict[repr(ghost_root_cell)] = ghost_root_cell.N

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

        RB, RBN, boundary_name_color_dict, pb_text = \
            cscg.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]
        reo_db = cscg.domain.regions.edges_on_domain_boundaries

        for rn in cscg.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = cscg.domain.regions.map[rn][ei]
                    if show_boundaries:
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                linewidth=3)
                    # noinspection PyUnresolvedReferences
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                            linewidth=0.75)

                if RBN[rn][ei] is None:
                    pass
                else:
                    if show_boundaries:
                        bn = cscg.domain.regions.map[rn][ei]
                        if bn in pb_text:
                            # noinspection PyUnresolvedReferences
                            ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                    '$<' + pb_text[bn] + '>$',
                                    c=boundary_name_color_dict[bn], ha='center', va='center')
                        else:
                            # noinspection PyUnresolvedReferences
                            ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                    '$<$' + bn + '$>$',
                                    c=boundary_name_color_dict[bn], ha='center', va='center')

        #---------------------- save to --------------------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')
        plt.close()






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
