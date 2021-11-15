# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Wed May 22 16:27:46 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
from root.config import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from SCREWS.frozen import FrozenOnly


class _2dCSCG_Domain_Visualize(FrozenOnly):
    """ """
    def __init__(self, domain):
        """ """
        assert domain.__class__.__name__ == '_2dCSCG_Domain', " <MeshMatchChecker> "
        assert domain.ndim == 2, " <MeshMatchChecker> "
        self._domain_ = domain
        self._freeze_self_()

    def __call__(self, **kwargs):
        if rAnk != mAster_rank: return
        return self.matplot(**kwargs)

    def matplot(self, show_region_boundary=True, usetex=True, corlormap='tab10',
        density=1000, xlim=None, ylim=None, show_region_names=True, show_boundary_names=True,
        labelsize=15, ticksize=15, fontsize=12, do_plot=True, saveto=None,
        domain_boundary_linewidth=3, region_linewidth=0.8, ):
        """

        :param show_region_boundary:
        :param usetex:
        :param corlormap:
        :param density:
        :param xlim:
        :param ylim:
        :param show_region_names:
        :param show_boundary_names:
        :param labelsize:
        :param ticksize:
        :param fontsize:
        :param do_plot:
        :param saveto:
        :param domain_boundary_linewidth:
        :param region_linewidth:
        :return:
        """
        if rAnk != mAster_rank: return
        density = int(np.ceil(density / self._domain_.regions.num / 4))
        if density > 100: density = 100
        if density < 30: density = 30
        o = np.linspace(0, 1, density)  # plot density
        O = np.zeros(density)
        I = np.ones(density)

        RB = {}  # region boundaries
        for rn in self._domain_.regions.names:
            RB[rn] = [None, None, None, None]
            for ei in range(4):
                if ei == 0:  # U
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(O, o)
                elif ei == 1:  # S
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(I, o)
                elif ei == 2:  # L
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(o, O)
                elif ei == 3:  # R
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(o, I)
                else:
                    raise Exception()
        if not do_plot:
            return RB

        if show_region_names:
            RC = {}
            RE = {}
            for rn in self._domain_.regions.names:
                RC[rn] = self._domain_.regions(rn).interpolation.mapping([0.5, ], [0.5, ])
                RE[rn] = [None, None, None, None]
                for ei in range(4):
                    if ei == 0:  # U
                        RE[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.15, ], [0.5, ])
                    elif ei == 1:  # D
                        RE[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.85, ], [0.5, ])
                    elif ei == 2:  # L
                        RE[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.5, ], [0.15, ])
                    elif ei == 3:  # R
                        RE[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.5, ], [0.85, ])
                    else:
                        raise Exception()

        reo_db = self._domain_.regions.edges_on_domain_boundaries

        if show_boundary_names:
            RBN = {}
            for rn in self._domain_.regions.names:
                RBN[rn] = [None, None, None, None]
                for ei in range(4):
                    if reo_db[rn][ei] == 1:
                        if ei == 0:  # U
                            RBN[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0, ], [0.5, ])
                        elif ei == 1:  # D
                            RBN[rn][ei] = self._domain_.regions(rn).interpolation.mapping([1, ], [0.5, ])
                        elif ei == 2:  # L
                            RBN[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.5, ], [0, ])
                        elif ei == 3:  # R
                            RBN[rn][ei] = self._domain_.regions(rn).interpolation.mapping([0.5, ], [1, ])
                        else:
                            raise Exception()

        boundaries_numb = self._domain_.boundaries.num
        boundaries_name = self._domain_.boundaries.names
        bounbary_name_color_dict = dict()
        if boundaries_numb > 10 and corlormap=='tab10': corlormap = 'viridis'
        color = cm.get_cmap(corlormap, boundaries_numb)
        colors = []
        for j in range(boundaries_numb):
            colors.append(color(j))
        for j, bn in enumerate(boundaries_name):
            bounbary_name_color_dict[bn] = colors[j]

        plt.rc('text', usetex=usetex)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)
        for rn in self._domain_.regions.names:
            if show_region_names:
                ax.text(RC[rn][0], RC[rn][1], rn.replace('_', '\_'), c='k', ha='center', va='center')
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = self._domain_.regions.map[rn][ei]
                    # noinspection PyUnresolvedReferences
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=bounbary_name_color_dict[bn],
                            linewidth=domain_boundary_linewidth)
                    # noinspection PyUnresolvedReferences
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                            linewidth=0.25*domain_boundary_linewidth)
                else:
                    if show_region_boundary:
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], '--', color='k',
                                linewidth=region_linewidth)
                if show_region_names:
                    # noinspection PyUnresolvedReferences
                    ax.text(RE[rn][ei][0], RE[rn][ei][1],
                            {0: 'U', 1: 'D', 2: 'L', 3: 'R'}[ei], c='r', ha='center', va='center')

                if show_boundary_names:
                    if RBN[rn][ei] is None:
                        pass
                    else:
                        bn = self._domain_.regions.map[rn][ei]
                        # noinspection PyUnresolvedReferences
                        ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                '$<$' + bn + '$>$', fontsize=fontsize,
                                c=bounbary_name_color_dict[bn], ha='center', va='center')

        plt.xlabel(r"$x$", fontsize=labelsize)
        plt.ylabel(r"$y$", fontsize=labelsize)
        plt.tick_params(axis='both', which='both', labelsize=ticksize)
        plt.tight_layout()

        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')

        plt.show()
        return fig