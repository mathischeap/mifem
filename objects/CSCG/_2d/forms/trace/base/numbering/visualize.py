# -*- coding: utf-8 -*-


from root.config.main import *
from components.freeze.main import FrozenOnly
import matplotlib.pyplot as plt
from matplotlib import cm


class _2dCSCG_Trace_Numbering_Visualize(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._mesh_ = tf.mesh
        self._trace_ = tf.mesh.trace
        self._default_ = 'matplot'
        self._freeze_self_()

    def __call__(self, **kwargs):
        return getattr(self, self._default_)(**kwargs)

    def matplot(self, **kwargs):
        return getattr(self, f'_matplot_{self._tf_.__class__.__name__}_numbering')(**kwargs)


    def _matplot_trace_mesh_BASE_(self, ax, region_boundary=True, density=10000, usetex=False,
                                  show_element_numbering=True,
                                  corlormap='tab10', boundary_name_fontsize=12, title=True,
                                  show_boundary_names=True, domain_boundary_linewidth=3,
                                  region_boundary_linewidth=1, element_linewidth=0.6,
                                  element_color='gray'):
        """

        :param ax:
        :param region_boundary:
        :param density:
        :param usetex:
        :param show_element_numbering:
        :param corlormap:
        :param boundary_name_fontsize:
        :param title:
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
        TEC_P = dict() # for finding  the position of the other mesh element edge of a periodic trace element.
        for i in self._trace_.elements:
            tei = self._trace_.elements[i]
            TED[i] = tei.coordinate_transformation.mapping(o)
            TEC[i] = tei.coordinate_transformation.mapping(c)
            if tei.IS_on_periodic_boundary:
                # for finding  the position of the other mesh element edge of a periodic trace element.
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
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)

            for i in ted:
                ax.plot(*ted[i], color=element_color, linewidth=element_linewidth)
                if show_element_numbering:
                    ax.text(*tec[i], "@${}$".format(i),
                             color = 'k', alpha=0.5,
                            ha='center', va='center')
                    if i in tec_p:
                        ax.text(*tec_p[i], "@${}$".format(i),
                             color = 'k', alpha=0.5,
                                ha='center', va='center')

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
                                        '$<' + pb_text[bn] + '>$', fontsize=boundary_name_fontsize,
                                        c=bounbary_name_color_dict[bn], ha='center', va='center')
                            else:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<$' + bn + '$>$', fontsize=boundary_name_fontsize,
                                        c=bounbary_name_color_dict[bn], ha='center', va='center')
            #______ show the element numbering ____________________________________________
            for rn in mesh.domain.regions.names:
                eccrn = element_center_coordinates[rn]
                AEGN = mesh.___PRIVATE_generate_ALL_element_global_numbering___()
                gnrn = AEGN[rn]
                for i in range(mesh.elements.layout[rn][0]):
                    for j in range(mesh.elements.layout[rn][1]):
                        ax.text(eccrn[0][i,j], eccrn[1][i,j], "$e{}$".format(gnrn[i,j]),
                                color = 'r',
                                ha='center', va='center')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$y$')
            if title is True:
                title = f"numbering of {self._tf_.orientation}-{self._tf_.k}-trace-form: " + \
                        f"'{self._tf_.standard_properties.name}'"
                plt.title(title)
            elif title is False:
                pass
            else:
                plt.title(title)

            # plt.tight_layout()
            return cm.get_cmap('Dark2', 8)




    def _matplot__1Trace_Outer_numbering(self, saveto=None, **kwargs):
        """"""
        assert self._tf_.k == 1
        nodes = self._tf_.space.nodes
        nx, ny = nodes
        cx = (nx[1:] + nx[:-1]) / 2
        cy = (ny[1:] + ny[:-1]) / 2

        elements = self._trace_.elements
        MAPPING = dict()
        GATHERING = dict()
        MP_P = dict() # for finding  the position of the other mesh element edge of a periodic trace element.
        _ = self._tf_.numbering.trace_element_wise
        for i in elements:
            ele = elements[i]
            if ele.CHARACTERISTIC_edge in 'UD':
                MAPPING[i] = ele.coordinate_transformation.mapping(cy)
            elif ele.CHARACTERISTIC_edge in 'LR':
                MAPPING[i] = ele.coordinate_transformation.mapping(cx)
            else:
                raise Exception()

            if ele.IS_on_periodic_boundary:
                # for finding  the position of the other mesh element edge of a periodic trace element.
                MP_P[i] = ele.NON_CHARACTERISTIC_position

            GATHERING[i] = self._tf_.numbering.trace_element_wise[i].full_vector

        MAPPING = COMM.gather(MAPPING, root=MASTER_RANK)
        GATHERING = COMM.gather(GATHERING, root=MASTER_RANK)
        MP_P = COMM.gather(MP_P, root=MASTER_RANK)


        if RANK == MASTER_RANK:
            mapping = dict()
            gathering = dict()
            for i, MPi in enumerate(MAPPING):
                mapping.update(MPi)
                gathering.update(GATHERING[i])

            mp_p = dict() # for finding  the position of the other mesh element edge of a periodic trace element.
            for mppi in MP_P:
                mp_p.update(mppi)
        else:
            mp_p = None

        mp_p = COMM.bcast(mp_p, root=MASTER_RANK)
        for pi in mp_p:
            position = mp_p[pi][-1]
            if position in 'UD':
                mp_p[pi] = self._trace_.elements.DO_compute_mapping_of_trace_at_position(mp_p[pi], cy)
            elif position in 'LR':
                mp_p[pi] = self._trace_.elements.DO_compute_mapping_of_trace_at_position(mp_p[pi], cx)
            else:
                raise Exception()


        if RANK == MASTER_RANK:
            fig, ax = plt.subplots(figsize=(15,9))
        else:
            ax = None

        LN_colors = self._matplot_trace_mesh_BASE_(ax, **kwargs)

        if RANK == MASTER_RANK:
            # .. now, we attach the numbering of 1-trace-form to the fig.
            for k in mapping: # go through all trace-elements (kth).
                mpk = mapping[k]
                gtk = gathering[k]
                ck = LN_colors(k % 8)
                for i, text in enumerate(gtk):
                    x = mpk[0][i]
                    y = mpk[1][i]
                    ax.text(x, y, str(text), c=ck, va='center', ha='center')

                    if k in mp_p:
                        x = mp_p[k][0][i]
                        y = mp_p[k][1][i]
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')
            #... SAVE TO ...
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
            #......
            plt.show()
            return fig



    def _matplot__1Trace_Inner_numbering(self, **kwargs):
        return self._matplot__1Trace_Outer_numbering(**kwargs)