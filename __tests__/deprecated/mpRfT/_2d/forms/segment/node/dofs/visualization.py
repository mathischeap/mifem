# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 6:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.freeze.base import FrozenOnly
import matplotlib.pyplot as plt
import matplotlib


class mpRfT2_NSgF_Dofs_Visualization(FrozenOnly):
    """"""

    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._t_ = dofs._t_
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    def matplot(self, density=None, color_space=True,
        saveto=None, usetex=False, labelsize=12, ticksize=12,
        ):
        """

        Parameters
        ----------
        density
        color_space
        saveto
        usetex
        labelsize
        ticksize

        Returns
        -------

        """
        mesh = self._t_.mesh
        CPD = dict()
        color_Space_dict = dict()
        for ci in mesh: # ao through all local cell indices.
            cell = mesh[ci]
            assert cell.___isroot___
            CPD[repr(cell)] = cell.coordinate_transformation.___PRIVATE_plot_data___(density=density)

            if color_space:
                color_Space_dict[repr(cell)] = cell.N

        CPD = cOmm.gather(CPD, root=mAster_rank)

        if color_space:
            color_Space_dict = cOmm.gather(color_Space_dict, root=mAster_rank)

        ntype = self._t_.ntype
        N = self._t_.N
        Nodes_dict = dict()
        for seg in mesh.segments:
            n = N[seg]

            if ntype == 'Gauss':
                nodes = mesh.space.Gauss(n)[0]

            else:
                raise NotImplementedError(f"not implemented for ntype={ntype}.")

            nodes = seg.coordinate_transformation.mapping(nodes)

            Nodes_dict[seg.__repr__()] = nodes

        Nodes_dict = cOmm.gather(Nodes_dict, root=mAster_rank)

        if rAnk != mAster_rank: return

        ___ = dict()
        for _ in CPD:
            ___.update(_)
        CPD = ___

        ___ = dict()
        for _ in Nodes_dict:
            ___.update(_)
        Nodes_dict = ___

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

        for srp in Nodes_dict:
            x, y = Nodes_dict[srp]
            plt.scatter(x, y, marker='s', s=4)

        #---------------------- save to ----------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')

        plt.close()






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/dofs/visualization.py
    pass
