# -*- coding: utf-8 -*-

from root.config import *
import matplotlib.pyplot as plt
from matplotlib import cm
from SCREWS.frozen import FrozenOnly


class _2dCSCG_Trace(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = _2dCSCG_Trace_Elements(self)
        self._visualize_ = _2dCSCG_Trace_Visualize(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.elements.RESET_cache()

    @property
    def elements(self):
        return self._elements_

    @property
    def visualize(self):
        return self._visualize_




class _2dCSCG_Trace_Elements(FrozenOnly):
    def __init__(self, trace):
        self._trace_ = trace
        self._mesh_ = trace._mesh_
        self.___PRIVATE_generating_trace_elements___()
        self.RESET_cache()
        self._freeze_self_()


    def RESET_cache(self):
        self._type_amount_dict_ = None

    def ___DO_find_type_and_amount_numbered_before___(self):
        """
        :return: A dictionary. For example, ``{..., 65: [32, 33,], ...}``, it means
            we have 32 'UD', 33 'LR' trace elements numbered before 65. We can see that
            32 + 33 = 65 (0, 1, 2, ..., 64).
        :rtype: dict
        """
        if self._type_amount_dict_ is None:
            if rAnk == 0:
                POOL = dict()
            else:
                POOL = cOmm.recv(source=rAnk - 1, tag=rAnk)

            type_amount_dict = dict()

            INDICES = list(self._elements_.keys())
            INDICES = sorted(INDICES)

            for i in INDICES:
                if i == 0: POOL[0] = np.array([0,0])

                tei = self[i]
                cs = tei.CHARACTERISTIC_edge
                POOL_i_p_1 = POOL[i].copy()
                if cs in 'UD':
                    POOL_i_p_1[0] += 1
                elif cs in 'LR':
                    POOL_i_p_1[1] += 1
                else:
                    raise Exception()
                if i+1 in POOL:
                    assert np.all(POOL_i_p_1 == POOL[i+1])
                else:
                    if i == self.GLOBAL_num - 1:
                        assert np.sum(POOL_i_p_1) == self.GLOBAL_num, "Something is wrong."
                    else:
                        POOL[i+1] = POOL_i_p_1

                type_amount_dict[i] = POOL[i]

                if not tei.IS_shared_by_cores:
                    del POOL[i]
                else:
                    CORE = tei.shared_with_core
                    if CORE < rAnk:
                        del POOL[i]

            if rAnk == sIze - 1:
                pass
            else:
                cOmm.send(POOL, dest=rAnk+1, tag=rAnk+1)

            self._type_amount_dict_ = type_amount_dict

        return self._type_amount_dict_


    @property
    def map(self):
        return self._MAP_

    @property
    def num(self):
        """
        (int)

        :return:
        """
        return len(self._elements_)

    def __getitem__(self, item):
        return self._elements_[item]

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i

    def __len__(self):
        return self.num


    @property
    def GLOBAL_num(self):
        """
        (int) The total number of trace elements.

        :return:
        """
        return self._GLOBAL_num_

    def ___PRIVATE_generating_trace_elements___(self):
        self._elements_ = dict()
        elements_map = self._mesh_.elements.map
        sideNames = 'UDLR'
        sidePairs = {'U':'D', 'D':'U', 'L':'R', 'R':'L'}
        upesp = self._mesh_.___useful_periodic_element_edge_pairs___
        self._MAP_ = dict()

        if rAnk == 0: # first core
            cn = 0
            POOL = dict()
        else: # intermediate cores
            POOL, cn = cOmm.recv(source=rAnk-1, tag=rAnk)

        LOCAL_POOL = dict()
        for i in elements_map:
            self._MAP_[i] = [None for _ in range(4)]
            for j in range(4):
                side_1 = sideNames[j]
                position_1 = str(i) + side_1
                what = elements_map[i][j]
                if isinstance(what, str): #on domain boundary
                    position_2 = what
                    self._elements_[cn] = _2dCSCG_Trace_Element(
                        self, cn, position_1, position_2, position_1, ondb=True, onpb=False)
                    self._MAP_[i][j] = cn
                    cn += 1
                else:
                    side_2 = sidePairs[sideNames[j]]
                    position_2 = str(what) + side_2
                    pair_1 = str(i) + '-' + side_1 + '|' + side_2 + '-' + str(what)
                    pair_2 = str(what) + '-' + side_2 + '|' + side_1 + '-' + str(i)
                    onpb = True if pair_1 in upesp or pair_2 in upesp else False
                    if what > i:
                        self._elements_[cn] = _2dCSCG_Trace_Element(
                            self, cn, position_1, position_2, position_1, ondb=False, onpb=onpb)
                        self._MAP_[i][j] = cn
                        POOL[position_1] = cn
                        POOL[position_2] = cn
                        cn += 1
                    elif what == i:
                        assert onpb
                        if side_1 in 'UL':
                            self._elements_[cn] = _2dCSCG_Trace_Element(
                                self, cn, position_1, position_2, position_1, ondb=False, onpb=True)
                            self._MAP_[i][j] = cn
                            LOCAL_POOL[position_2] = cn
                            cn += 1
                        else:
                            self._MAP_[i][j] = LOCAL_POOL[position_1]
                            del LOCAL_POOL[position_1]
                    else:
                        assert position_1 in POOL and position_2 in POOL
                        assert POOL[position_1] == POOL[position_2]
                        tn = POOL[position_1]
                        self._MAP_[i][j] = tn
                        if tn in self._elements_:
                            pass
                        else:
                            self._elements_[tn] = _2dCSCG_Trace_Element(
                                self, tn, position_2, position_1, position_1, ondb=False, onpb=onpb)

                        del POOL[position_1]
                        del POOL[position_2]

        if rAnk == sIze - 1:
            assert len(POOL) == 0
        else:
            cOmm.send([POOL, cn], dest=rAnk+1, tag=rAnk+1)

        cOmm.barrier()

        cn = cOmm.bcast(cn, root=sIze-1)
        self._GLOBAL_num_ = cn

        if saFe_mode:
            for i in elements_map:
                for j in range(4):
                    assert self._MAP_[i][j] in self, f"Miss local trace element {self._MAP_[i][j]}"


    def DO_compute_mapping_of_trace_at_position(self, position, ep1d):
        """

        :param position: Like "10-U", then we compute for the trace element which is the upper edge of
            10th mesh element.
        :param ep1d: Just like the ep1d in element.ct.
        :return: Only return to master core.
        """
        edge = position[-1]
        element = int(position[:-1])

        if element in self._mesh_.elements:
            ep = self.___generate_full_ep___(ep1d, edge)
            xy = self._mesh_.elements[element].coordinate_transformation.mapping(*ep)
        else:
            xy = None

        xy = cOmm.gather(xy, root=mAster_rank)
        if rAnk == mAster_rank:
            assert xy.count(None) == sIze - 1
            return next(item for item in xy if item is not None) # find the first item which is not None.

    @staticmethod
    def ___generate_full_ep___(ep1d, element_edge):
        """
        Given the 1d evaluation points, we generate the 2d coordinates according to
        the element_edge.

        Parameters
        ----------
        ep1d:
            The 1d evaluation points.
        element_edge:
            The mesh-element edge: U, D, L, or R.

        """
        assert np.ndim(ep1d) == 1, " <TraceElementCoordinateTransformation> "
        assert np.min(ep1d) >= -1 and np.max(ep1d) <= 1, " <TraceElementCoordinateTransformation> "
        if np.size(ep1d) > 1:
            assert all(np.diff(ep1d) > 0), " <TraceElementCoordinateTransformation> "
        shape = len(ep1d)
        if element_edge == 'U':
            ep = (-np.ones(shape), ep1d)
        elif element_edge == 'D':
            ep = (np.ones(shape), ep1d)
        elif element_edge == 'L':
            ep = (ep1d, -np.ones(shape))
        elif element_edge == 'R':
            ep = (ep1d, np.ones(shape))
        else:
            raise Exception()
        return ep


class _2dCSCG_Trace_Element(FrozenOnly):
    """

    :param trace_elements:
    :param i:
    :param position_1:
    :param position_2:
    :param cp: characteristic position
    :param ondb:
    :param onpb:
    """
    def __init__(self, trace_elements, i, position_1, position_2, cp, ondb=False, onpb=False):
        self._elements_ = trace_elements
        self._mesh_ = trace_elements._mesh_
        self._i_ = i
        self._p1_ = position_1
        self._p2_ = position_2
        self._cp_ = cp
        if position_1 == cp:
            self._ncp_ = position_2
        elif position_2 == cp:
            self._ncp_ = position_1
        else:
            raise Exception()
        self._ondb_ = ondb
        self._onpb_ = onpb
        assert self.CHARACTERISTIC_element in self._elements_._mesh_.elements, \
            "CHARACTERISTIC_element must be int the same core."
        if self.IS_on_mesh_boundary:
            assert self.NON_CHARACTERISTIC_position[0] not in '1234567890'
        self._ct_ = _2dCSCG_Trace_Element_CoordinateTransformation(self)
        self._freeze_self_()

    @property
    def positions(self):
        return self._p1_, self._p2_

    @property
    def coordinate_transformation(self):
        return self._ct_

    @property
    def normal_direction(self):
        """"""
        if self._p1_[-1] in 'NS':
            return 'NS'
        elif self._p1_[-1] in 'WE':
            return 'WE'
        elif self._p1_[-1] in 'BF':
            return 'BF'
        else:
            raise Exception()

    @property
    def NON_CHARACTERISTIC_position(self):
        return self._ncp_

    @property
    def CHARACTERISTIC_position(self):
        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        return int(self._cp_[:-1])
    @property
    def CHARACTERISTIC_edge(self):
        return self._cp_[-1]

    @property
    def i(self):
        return self._i_

    @property
    def IS_on_mesh_boundary(self):
        return self._ondb_

    @property
    def IS_on_periodic_boundary(self):
        return self._onpb_

    @property
    def IS_shared_by_cores(self):
        """"""
        if self.IS_on_mesh_boundary:
            return False
        else:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements and \
                int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                return False
            else:
                return True

    @property
    def shared_with_core(self):
        if self.IS_shared_by_cores:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.DO.FIND_slave_of_element(int(self._p2_[:-1]))
            elif int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.DO.FIND_slave_of_element(int(self._p1_[:-1]))
            else:
                raise Exception()
            assert CORE != rAnk
            return CORE
        else:
            return None



class _2dCSCG_Trace_Element_CoordinateTransformation(FrozenOnly):
    def __init__(self, te):
        self._te_ = te
        self._freeze_self_()


    def mapping(self, evaluation_points):
        """
        The local mapping.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_edge
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side)
        x, y = self._te_._mesh_.elements[i].coordinate_transformation.mapping(*ep)
        return x, y

    def Jacobian_matrix(self, evaluation_points):
        """
        The local Jacobian matrix.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_edge = self._te_.CHARACTERISTIC_edge
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_edge)
        J = self._te_._mesh_.elements[i].coordinate_transformation.Jacobian_matrix(*ep)
        if element_edge in 'UD':
            return J[0][1], J[1][1]
        elif element_edge in 'LR':
            return J[0][0], J[1][0]
        else:
            raise Exception()

    def metric_matrix(self, ep1d):
        """ The entries of metric_matrix is normally denoted as g_{i,j}. """
        J = self.Jacobian_matrix(ep1d)
        Gi = J[0]**2 + J[1]**2
        return Gi

    def metric(self, ep1d):
        """ g, which should be det(metric_matrix): det(G). """
        return self.metric_matrix(ep1d)






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
        domain_boundary_linewidth=3, region_boundary_linewidth=1, element_linwidth=0.4,
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
        :param element_linwidth:
        :param element_color:
        :return:
        """
        mesh = self._trace_._mesh_
        density = int(np.ceil(density/self._trace_.elements.GLOBAL_num))
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
            if tei.IS_on_periodic_boundary:
                TEC_P[i] = tei.NON_CHARACTERISTIC_position

        TED = cOmm.gather(TED, root=mAster_rank)
        TEC = cOmm.gather(TEC, root=mAster_rank)
        TEC_P = cOmm.gather(TEC_P, root=mAster_rank)


        if rAnk == mAster_rank:
            # for finding  the position of the other mesh element edge of a periodic trace element.
            tec_p = dict()
            for TEcp_i in TEC_P:  tec_p.update(TEcp_i)
        else:
            tec_p = None
        tec_p = cOmm.bcast(tec_p, root=mAster_rank)
        for pi in tec_p:
            tec_p[pi] = self._trace_.elements.DO_compute_mapping_of_trace_at_position(tec_p[pi], c)

        if rAnk == mAster_rank:
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
                ax.plot(*ted[i], color=element_color, linewidth=element_linwidth)
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

