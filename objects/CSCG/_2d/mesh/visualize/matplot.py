# -*- coding: utf-8 -*-
from root.config.main import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm

from components.freeze.main import FrozenOnly


class _2dCSCG_Mesh_Visualize_Matplot(FrozenOnly):
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', " <MeshVisualize> "
        assert mesh.ndim == 2, " <MeshVisualize> "
        self._mesh_ = mesh
        self._domain_ = mesh.domain
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self._matplot_mesh_(*args, **kwargs)

    def element_division(
            self, density=50000, usetex=False, saveto=None,
            xlim=None, ylim=None, labelsize=15, ticksize=15, element_linewidth=0.4,
            element_color='red', corlormap=None
    ):
        """plot element division."""
        if RANK != MASTER_RANK:
            return

        density = int(np.ceil(density / self._mesh_.elements.global_num))
        max_element_layout = 0
        for rn in self._mesh_.domain.regions.names:
            if np.max(self._mesh_.elements.layout[rn]) > max_element_layout:
                max_element_layout = np.max(self._mesh_.elements.layout[rn])
        if density > 30 * max_element_layout:
            density = 30 * max_element_layout
        if density < 3 * max_element_layout:
            density = 3 * max_element_layout

        o = np.linspace(0, 1, density)  # plot density
        _I = np.ones(density)
        # ______________________________________ line data _____________________________
        RI = {}  # data for regions internal lines
        for rn in self._mesh_.domain.regions.names:
            RI[rn] = ([], [])  # ([dy_lines], [dx_lines])
            # ____ compute dy lines ___________________________________________
            for x in self._mesh_.elements.spacing[rn][0][:]:
                RI[rn][0].append(self._mesh_.domain.regions(rn).interpolation.mapping(x * _I, o))
            # ____ compute dx lines ___________________________________________
            for y in self._mesh_.elements.spacing[rn][1][:]:
                RI[rn][1].append(self._mesh_.domain.regions(rn).interpolation.mapping(o, y * _I))
            # --------------------------------------------------------------------------

        # _____________ text: element numbering data ___________________________________
        element_center_coordinates = {}
        for rn in self._mesh_.domain.regions.names:
            element_center_coordinate_xi = (self._mesh_.elements.spacing[rn][0][:-1]
                                            + self._mesh_.elements.spacing[rn][0][1:]) / 2
            element_center_coordinate_eta = (self._mesh_.elements.spacing[rn][1][:-1]
                                             + self._mesh_.elements.spacing[rn][1][1:]) / 2
            element_center_coordinate_eta, element_center_coordinate_xi = \
                np.meshgrid(element_center_coordinate_eta, element_center_coordinate_xi)
            element_center_coordinates[rn] = \
                self._mesh_.domain.regions(rn).interpolation.mapping(
                        element_center_coordinate_xi,
                        element_center_coordinate_eta)

        # ___________ do the plot ______________________________________________________
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
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        for rn in self._mesh_.domain.regions.names:
            for dy_lines in RI[rn][0]:  # plot dy mesh lines
                ax.plot(dy_lines[0], dy_lines[1], color=element_color, linewidth=element_linewidth)
            for dx_lines in RI[rn][1]:  # plot dx mesh lines
                ax.plot(dx_lines[0], dx_lines[1], color=element_color, linewidth=element_linewidth)

        RB = self._mesh_.domain.visualize.matplot(
                density=4*150*self._mesh_.domain.regions.num, do_plot=False)
        reodb = self._mesh_.domain.regions.edges_on_domain_boundaries
        for rn in self._mesh_.domain.regions.names:
            for ei in range(4):
                if reodb[rn][ei] == 1:  # plot the domain boundary
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k', linewidth=element_linewidth*3)
                else:
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='b', linewidth=element_linewidth*2)

        # prepare colors for different cores...
        if corlormap is None:
            corlormap = 'tab20'
        color = cm.get_cmap(corlormap, SIZE)
        COLOR = dict()
        for i in range(SIZE):
            COLOR[i] = color(i)

        for rn in self._mesh_.domain.regions.names:
            eccrn = element_center_coordinates[rn]
            AEGN = self._mesh_.___PRIVATE_generate_ALL_element_global_numbering___()
            gnrn = AEGN[rn]
            for i in range(self._mesh_.elements.layout[rn][0]):
                for j in range(self._mesh_.elements.layout[rn][1]):

                    C = self._mesh_.do.find.slave_of_element(gnrn[i, j])

                    if C == MASTER_RANK:  # for the master core, we box the element numbering.
                        plt.text(eccrn[0][i, j], eccrn[1][i, j], "${}$".format(gnrn[i, j]),
                                 bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 0},
                                 color=COLOR[C], fontsize=11, ha='center', va='center')
                    else:
                        plt.text(eccrn[0][i, j], eccrn[1][i, j], "${}$".format(gnrn[i, j]),
                                 color=COLOR[C], fontsize=11, ha='center', va='center')

        plt.tight_layout()
        # __________ SAVE TO ___________________________________________________________
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
        else:
            plt.show()

        plt.close('all')
        # -------------------------------------------------------------------------------
        return fig

    def ___PRIVATE_DO_generate_boundary_data___(self, density, usetex=True):
        """ To use this, do:

        RB, RBN, boundary_name_color_dict, pb_text = \
            self._mesh_.visualize.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]

        reo_db = self._mesh_.domain.regions.edges_on_domain_boundaries

        for rn in self._mesh_.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = self._mesh_.domain.regions.map[rn][ei]
                    # noinspection PyUnresolvedReferences
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                            linewidth=domain_boundary_linewidth)
                    # noinspection PyUnresolvedReferences
                    ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                            linewidth=0.25*domain_boundary_linewidth)

                if RBN[rn][ei] is None:
                    pass
                else:
                    bn = self._mesh_.domain.regions.map[rn][ei]
                    if bn in pb_text:
                        # noinspection PyUnresolvedReferences
                        ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                '$<' + pb_text[bn] + '>$', fontsize=boundary_name_fontsize,
                                c=boundary_name_color_dict[bn])
                    else:
                        # noinspection PyUnresolvedReferences
                        ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                '$<$' + bn + '$>$', fontsize=boundary_name_fontsize,
                                c=boundary_name_color_dict[bn])

        """

        o = np.linspace(0, 1, density)  # plot density
        _O = np.zeros(density)
        _I = np.ones(density)
        RB = {}  # regions boundaries
        for rn in self._domain_.regions.names:
            RB[rn] = [None, None, None, None]
            for ei in range(4):
                if ei == 0:  # U
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(_O, o)
                elif ei == 1:  # S
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(_I, o)
                elif ei == 2:  # L
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(o, _O)
                elif ei == 3:  # R
                    RB[rn][ei] = self._domain_.regions(rn).interpolation.mapping(o, _I)
                else:
                    raise Exception()
        RBN = {}
        reo_db = self._domain_.regions.edges_on_domain_boundaries
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
        boundary_name_color_dict: dict = dict()
        if boundaries_numb > 10:
            corlormap = 'viridis'
        else:
            corlormap = 'tab10'
        color = cm.get_cmap(corlormap, boundaries_numb)
        colors = []
        for j in range(boundaries_numb):
            colors.append(color(j))
        for j, bn in enumerate(boundaries_name):
            boundary_name_color_dict[bn] = colors[j]

        # .. take care periodic boundaries
        DI = self._mesh_.domain.domain_input
        pbp = DI.periodic_boundary_pairs

        pbs = DI.periodic_boundaries
        pb_text = dict()
        if pbp == dict():  # no periodic boundaries, lets just pass.
            assert pbs == set()
        else:
            for pair in pbp:
                pb1, pb2 = pair.split('=')
                ptype = pbp[pair]
                boundary_name_color_dict[pb2] = boundary_name_color_dict[pb1]

                if usetex:
                    pb_text[pb1] = r'\mathrm{%s}' % pb1 + \
                                   r'\stackrel{\mathrm{%s}}{=}' % ptype + r'\mathrm{%s}' % pb2
                    pb_text[pb2] = r'\mathrm{%s}' % pb2 + \
                                   r'\stackrel{\mathrm{%s}}{=}' % ptype + r'\mathrm{%s}' % pb1
                else:
                    pb_text[pb1] = r'\mathrm{%s}' % pb1 + \
                                   r'\genfrac{}{}{0}{}{%s}{=}' % ptype + r'\mathrm{%s}' % pb2
                    pb_text[pb2] = r'\mathrm{%s}' % pb2 + \
                                   r'\genfrac{}{}{0}{}{%s}{=}' % ptype + r'\mathrm{%s}' % pb1

        return RB, RBN, boundary_name_color_dict, pb_text

    def _matplot_mesh_(
            self, region_boundary=True, density=6000, usetex=False,
            show_numbering=False,
            saveto=None, pad_inches=0,
            corlormap='tab10', fontsize=12,
            xlim=None, ylim=None, xticks=None, yticks=None,
            labelsize=15, ticksize=15, show_boundary_names=True,
            highlight_domain_boundary=True,
            domain_boundary_linewidth=3, region_boundary_linewidth=0.8, element_linewidth=0.4,
            element_color='red', top_spine=True, bottom_spine=True, left_spine=True, right_spine=True,
    ):
        """

        Parameters
        ----------
        region_boundary
        density
        usetex
        show_numbering
        saveto
        pad_inches
        corlormap
        fontsize
        xlim
        ylim
        xticks
        yticks
        labelsize
        ticksize
        show_boundary_names
        highlight_domain_boundary
        domain_boundary_linewidth
        region_boundary_linewidth
        element_linewidth
        element_color
        top_spine
        bottom_spine
        left_spine
        right_spine

        Returns
        -------

        """
        if RANK != MASTER_RANK:
            return

        density = int(np.ceil(density / self._mesh_.elements.global_num))
        max_element_layout = 0
        for rn in self._mesh_.domain.regions.names:
            if np.max(self._mesh_.elements.layout[rn]) > max_element_layout:
                max_element_layout = np.max(self._mesh_.elements.layout[rn])
        if density > 30 * max_element_layout:
            density = 30 * max_element_layout
        if density < 3 * max_element_layout:
            density = 3 * max_element_layout
        RB = self._mesh_.domain.visualize.matplot(
                density=4*150*self._mesh_.domain.regions.num, do_plot=False)
        o = np.linspace(0, 1, density)  # plot density
        _I = np.ones(density)
        # ______________________________________ line data _____________________________
        RI = {}  # data for regions internal lines
        for rn in self._mesh_.domain.regions.names:
            RI[rn] = ([], [])   # ([dy_lines], [dx_lines])
            # ____ compute dy lines ___________________________________________
            for x in self._mesh_.elements.spacing[rn][0][:]:
                RI[rn][0].append(self._mesh_.domain.regions(rn).interpolation.mapping(x * _I, o))
            # ____ compute dx lines ___________________________________________
            for y in self._mesh_.elements.spacing[rn][1][:]:
                RI[rn][1].append(self._mesh_.domain.regions(rn).interpolation.mapping(o, y * _I))
            # --------------------------------------------------------------------------
        # _____________ text: element numbering data ___________________________________
        if show_numbering:
            element_center_coordinates = {}
            for rn in self._mesh_.domain.regions.names:
                element_center_coordinate_xi = (self._mesh_.elements.spacing[rn][0][:-1]
                                                + self._mesh_.elements.spacing[rn][0][1:]) / 2
                element_center_coordinate_eta = (self._mesh_.elements.spacing[rn][1][:-1]
                                                 + self._mesh_.elements.spacing[rn][1][1:]) / 2
                element_center_coordinate_eta, element_center_coordinate_xi = \
                    np.meshgrid(element_center_coordinate_eta, element_center_coordinate_xi)
                element_center_coordinates[rn] = \
                    self._mesh_.domain.regions(rn).interpolation.mapping(
                            element_center_coordinate_xi,
                            element_center_coordinate_eta)

        reodb = self._mesh_.domain.regions.edges_on_domain_boundaries
        # _____________ text: element numbering data ___________________________________
        if show_boundary_names:
            RBN = {}
            for rn in self._mesh_.domain.regions.names:
                RBN[rn] = [None, None, None, None]
                for ei in range(4):
                    if reodb[rn][ei] == 1:
                        if ei == 0:  # U
                            RBN[rn][ei] = self._mesh_.domain.regions(rn).interpolation.mapping([0, ], [0.5, ])
                        elif ei == 1:  # D
                            RBN[rn][ei] = self._mesh_.domain.regions(rn).interpolation.mapping([1, ], [0.5, ])
                        elif ei == 2:  # L
                            RBN[rn][ei] = self._mesh_.domain.regions(rn).interpolation.mapping([0.5, ], [0, ])
                        elif ei == 3:  # R
                            RBN[rn][ei] = self._mesh_.domain.regions(rn).interpolation.mapping([0.5, ], [1, ])
                        else:
                            raise Exception()

        # _____ get personal color for boundaries ________________________________________
        boundaries_numb = self._mesh_.domain.boundaries.num
        boundaries_name = self._mesh_.domain.boundaries.names
        boundary_name_color_dict = dict()
        if boundaries_numb > 10 and corlormap == 'tab10':
            corlormap = 'viridis'
        color = cm.get_cmap(corlormap, boundaries_numb)
        colors = []
        for j in range(boundaries_numb):
            colors.append(color(j))
        for j, bn in enumerate(boundaries_name):
            boundary_name_color_dict[bn] = colors[j]

        # .. take care periodic boundaries
        DI = self._mesh_.domain.domain_input
        pbp = DI.periodic_boundary_pairs
        pbs = DI.periodic_boundaries
        pb_text = dict()
        if pbp == dict():  # no periodic boundaries, lets just pass.
            assert pbs == set()
        else:
            for pair in pbp:
                pb1, pb2 = pair.split('=')
                ptype = pbp[pair]
                # noinspection PyUnresolvedReferences
                boundary_name_color_dict[pb2] = boundary_name_color_dict[pb1]

                if usetex:
                    pb_text[pb1] = r'\mathrm{%s}' % pb1 + \
                                   r'\stackrel{\mathrm{%s}}{=}' % ptype + r'\mathrm{%s}' % pb2
                    pb_text[pb2] = r'\mathrm{%s}' % pb2 + \
                                   r'\stackrel{\mathrm{%s}}{=}' % ptype + r'\mathrm{%s}' % pb1
                else:
                    pb_text[pb1] = r'\mathrm{%s}' % pb1 + \
                                   r'\genfrac{}{}{0}{}{%s}{=}' % ptype + r'\mathrm{%s}' % pb2
                    pb_text[pb2] = r'\mathrm{%s}' % pb2 + \
                                   r'\genfrac{}{}{0}{}{%s}{=}' % ptype + r'\mathrm{%s}' % pb1

        # ___________ do the plot ______________________________________________________
        if saveto is not None:
            matplotlib.use('Agg')
        plt.rcParams.update({
            "text.usetex": usetex,
            "font.family": "Times New Roman"
        })
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(top_spine)
        ax.spines['right'].set_visible(right_spine)
        ax.spines['left'].set_visible(left_spine)
        ax.spines['bottom'].set_visible(bottom_spine)
        if labelsize == 0:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            plt.xlabel(r"$x$", fontsize=labelsize)
            plt.ylabel(r"$y$", fontsize=labelsize)
        plt.tick_params(axis='both', which='both', labelsize=ticksize)

        if xticks is False:
            plt.xticks([])
        elif xticks is not None:
            plt.xticks(xticks)

        if yticks is False:
            plt.yticks([])
        elif yticks is not None:
            plt.yticks(yticks)

        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        for rn in self._mesh_.domain.regions.names:
            for dy_lines in RI[rn][0]:  # plot dy mesh lines
                ax.plot(dy_lines[0], dy_lines[1], color=element_color, linewidth=element_linewidth)
            for dx_lines in RI[rn][1]:  # plot dx mesh lines
                ax.plot(dx_lines[0], dx_lines[1], color=element_color, linewidth=element_linewidth)
        for rn in self._mesh_.domain.regions.names:
            for ei in range(4):
                if reodb[rn][ei] == 1:  # plot the domain boundary
                    bn = self._mesh_.domain.regions.map[rn][ei]
                    if highlight_domain_boundary:
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color=boundary_name_color_dict[bn], linewidth=domain_boundary_linewidth)
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color='k', linewidth=0.1*domain_boundary_linewidth)
                    else:
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color='k', linewidth=2*element_linewidth)

                    # Not an error, we just plot a thin line in line
                else:
                    if region_boundary:  # plot the regions boundary
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='b', linewidth=region_boundary_linewidth)

                if show_boundary_names:
                    # noinspection PyUnboundLocalVariable
                    if RBN[rn][ei] is None:
                        pass
                    else:
                        bn = self._mesh_.domain.regions.map[rn][ei]
                        if bn in pb_text:
                            # noinspection PyUnresolvedReferences
                            ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                    '$<' + pb_text[bn] + '>$', fontsize=fontsize,
                                    c=boundary_name_color_dict[bn], ha='center', va='center')
                        else:
                            # noinspection PyUnresolvedReferences
                            ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                    '$<$' + bn + '$>$', fontsize=fontsize,
                                    c=boundary_name_color_dict[bn], ha='center', va='center')
        # ______ show the element numbering ____________________________________________
        if show_numbering:
            for rn in self._mesh_.domain.regions.names:
                # noinspection PyUnboundLocalVariable
                eccrn = element_center_coordinates[rn]
                AEGN = self._mesh_.___PRIVATE_generate_ALL_element_global_numbering___()
                gnrn = AEGN[rn]
                for i in range(self._mesh_.elements.layout[rn][0]):
                    for j in range(self._mesh_.elements.layout[rn][1]):
                        plt.text(eccrn[0][i, j], eccrn[1][i, j], "${}$".format(gnrn[i, j]),
                                 color='gray', fontsize=11, ha='center', va='center')

        plt.tight_layout()
        # __________ SAVE TO ___________________________________________________________
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight', pad_inches=pad_inches)
        else:
            plt.show()

        plt.close()
        # ------------------------------------------------------------------------------
        return fig
