# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import cm
from components.freeze.main import FrozenOnly
from root.config.main import *






class _2dCSCG_Trace_Visualize(FrozenOnly):
    def __init__(self, trace):
        self._trace_ = trace
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    def matplot(self, region_boundary=True, density=10000, usetex=False,
                show_element_numbering=True, element_numbering_fontsize=12,
                saveto=None, corlormap='tab10', fontsize=12,
                xlim=None, ylim=None, labelsize=15, ticksize=15,
                show_boundary_names=True,
                domain_boundary_linewidth=3, region_boundary_linewidth=1, element_linewidth=0.4,
                element_color='red'):
        """

        :param region_boundary:
        :param density:
        :param usetex:
        :param show_element_numbering:
        :param element_numbering_fontsize:
        :param saveto:
        :param corlormap:
        :param fontsize:
        :param xlim:
        :param ylim:
        :param labelsize:
        :param ticksize:
        :param show_boundary_names:
        :param domain_boundary_linewidth:
        :param region_boundary_linewidth:
        :param element_linewidth:
        :param element_color:
        :return:
        """
        mesh = self._trace_._mesh_
        density = int(np.ceil(density / self._trace_.elements.global_num))
        if density > 100: density = 100
        if density < 10: density = 10

        o = np.linspace(-1, 1, density)  # plot density
        c = np.array([0,])
        TED = dict()
        TEC = dict()
        TEC_P = dict()
        for i in self._trace_.elements:
            tei = self._trace_.elements[i]
            TED[i] = tei.coordinate_transformation.mapping(o)
            TEC[i] = tei.coordinate_transformation.mapping(c)
            if tei.whether.on_periodic_boundary:
                TEC_P[i] = tei.NON_CHARACTERISTIC_position

        TED = COMM.gather(TED, root=MASTER_RANK)
        TEC = COMM.gather(TEC, root=MASTER_RANK)
        TEC_P = COMM.gather(TEC_P, root=MASTER_RANK)


        if RANK == MASTER_RANK:
            # for finding  the position of the other mesh element edge of a periodic trace element.
            tec_p = dict()
            for TEcp_i in TEC_P:  tec_p.update(TEcp_i)
        else:
            tec_p = None
        tec_p = COMM.bcast(tec_p, root=MASTER_RANK)
        for pi in tec_p:
            tec_p[pi] = self._trace_.elements.DO_compute_mapping_of_trace_at_position(tec_p[pi], c)

        if RANK == MASTER_RANK:
            ted, tec = dict(), dict()
            for TEDi in TED:  ted.update(TEDi)
            for TECi in TEC:  tec.update(TECi)
            del TED, TEC

            RB = mesh.domain.visualize.matplot(
                 density=4*150*mesh.domain.regions.num, do_plot=False)
            #_____________ text: element numbering data ___________________________________
            if show_element_numbering:
                element_center_coordinates = {}
                for rn in mesh.domain.regions.names:
                    element_center_coordinate_xi = (mesh.elements.spacing[rn][0][:-1]
                                                  +mesh.elements.spacing[rn][0][1:]) / 2
                    element_center_coordinate_eta = (mesh.elements.spacing[rn][1][:-1]
                                                  +mesh.elements.spacing[rn][1][1:]) / 2
                    element_center_coordinate_eta, element_center_coordinate_xi = \
                        np.meshgrid(element_center_coordinate_eta, element_center_coordinate_xi)
                    element_center_coordinates[rn] = \
                        mesh.domain.regions(rn).interpolation.mapping(
                                element_center_coordinate_xi,
                                element_center_coordinate_eta)

            reodb = mesh.domain.regions.edges_on_domain_boundaries
            # _____________ text: element numbering data ___________________________________
            if show_boundary_names:
                RBN = {}
                for rn in mesh.domain.regions.names:
                    RBN[rn] = [None, None, None, None]
                    for ei in range(4):
                        if reodb[rn][ei] == 1:
                            if ei == 0:  # U
                                RBN[rn][ei] = mesh.domain.regions(rn).interpolation.mapping([0, ], [0.5, ])
                            elif ei == 1:  # D
                                RBN[rn][ei] = mesh.domain.regions(rn).interpolation.mapping([1, ], [0.5, ])
                            elif ei == 2:  # L
                                RBN[rn][ei] = mesh.domain.regions(rn).interpolation.mapping([0.5, ], [0, ])
                            elif ei == 3:  # R
                                RBN[rn][ei] = mesh.domain.regions(rn).interpolation.mapping([0.5, ], [1, ])
                            else:
                                raise Exception()

            #_____ get personal color for boundaries ________________________________________
            boundaries_numb = mesh.domain.boundaries.num
            boundaries_name = mesh.domain.boundaries.names
            bounbary_name_color_dict = dict()
            if boundaries_numb > 10 and corlormap=='tab10': corlormap = 'viridis'
            color = cm.get_cmap(corlormap, boundaries_numb)
            colors = []
            for j in range(boundaries_numb):
                colors.append(color(j))
            for j, bn in enumerate(boundaries_name):
                bounbary_name_color_dict[bn] = colors[j]

            # .. take care periodic boundaries ...
            DI = mesh.domain.domain_input
            pbp = DI.periodic_boundary_pairs
            pbs = DI.periodic_boundaries
            pb_text = dict()
            if pbp == dict(): # no periodic boundaries, lets just pass.
                assert pbs == set()
            else:
                for pair in pbp:
                    pb1, pb2 = pair.split('=')
                    ptype = pbp[pair]
                    bounbary_name_color_dict[pb2] = bounbary_name_color_dict[pb1]
                    if usetex:
                        pb_text[pb1] = '\mathrm{%s}'%pb1 + \
                                       '\stackrel{\mathrm{%s}}{=}'%ptype + '\mathrm{%s}'%pb2
                        pb_text[pb2] = '\mathrm{%s}'%pb2 + \
                                       '\stackrel{\mathrm{%s}}{=}'%ptype + '\mathrm{%s}'%pb1
                    else:
                        pb_text[pb1] = '\mathrm{%s}'%pb1 + \
                                       '\genfrac{}{}{0}{}{%s}{=}'%ptype + '\mathrm{%s}'%pb2
                        pb_text[pb2] = '\mathrm{%s}'%pb2 + \
                                       '\genfrac{}{}{0}{}{%s}{=}'%ptype + '\mathrm{%s}'%pb1

            #___________ do the plot ______________________________________________________
            plt.rc('text', usetex=usetex)
            fig, ax = plt.subplots(figsize=(15, 9))
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            plt.xlabel(r"$x$", fontsize=labelsize)
            plt.ylabel(r"$y$", fontsize=labelsize)
            plt.tick_params(axis='both', which='both', labelsize=ticksize)
            if xlim is not None: plt.xlim(xlim)
            if ylim is not None: plt.ylim(ylim)

            for i in ted:
                ax.plot(*ted[i], color=element_color, linewidth=element_linewidth)
                ax.text(*tec[i], "${}$".format(i),
                         color = 'k', fontsize=element_numbering_fontsize, ha='center', va='center')
                if i in tec_p:
                    ax.text(*tec_p[i], "${}$".format(i),
                            color='k', fontsize=element_numbering_fontsize, ha='center', va='center')

            for rn in mesh.domain.regions.names:
                for ei in range(4):
                    if reodb[rn][ei] == 1: # plot the domain boundary
                        bn = mesh.domain.regions.map[rn][ei]
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color=bounbary_name_color_dict[bn], linewidth=domain_boundary_linewidth)
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color='k', linewidth=0.1*domain_boundary_linewidth)
                    else:
                        if region_boundary: # plot the regions boundary
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='b',
                                    linewidth=region_boundary_linewidth)

                    if show_boundary_names:
                        if RBN[rn][ei] is None:
                            pass
                        else:
                            bn = mesh.domain.regions.map[rn][ei]
                            if bn in pb_text:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<' + pb_text[bn] + '>$', fontsize=fontsize,
                                        c=bounbary_name_color_dict[bn], ha='center', va='center')
                            else:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<$' + bn + '$>$', fontsize=fontsize,
                                        c=bounbary_name_color_dict[bn], ha='center', va='center')
            #______ show the element numbering ____________________________________________
            if show_element_numbering:
                for rn in mesh.domain.regions.names:
                    eccrn = element_center_coordinates[rn]
                    AEGN = mesh.___PRIVATE_generate_ALL_element_global_numbering___()
                    gnrn = AEGN[rn]
                    for i in range(mesh.elements.layout[rn][0]):
                        for j in range(mesh.elements.layout[rn][1]):
                            ax.text(eccrn[0][i,j], eccrn[1][i,j], "$e{}$".format(gnrn[i,j]),
                                    color = 'r', fontsize=element_numbering_fontsize,
                                    ha='center', va='center')

            plt.tight_layout()
            #__________ SAVE TO ___________________________________________________________
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
            #------------------------------------------------------------------------------
            plt.show()
            return fig

