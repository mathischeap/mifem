# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config.main import *
from components.freeze.main import FrozenOnly
import matplotlib.pyplot as plt
from matplotlib import cm


class _2dCSCG_Numbering_Visualize(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self._mesh_ = f.mesh
        self._default_ = 'matplot'
        self._freeze_self_()

    def __call__(self, **kwargs):
        return getattr(self, self._default_)(**kwargs)

    def matplot(self, **kwargs):
        return getattr(self, f'_matplot_{self._f_.__class__.__name__}_numbering')(**kwargs)

    def _matplot_mesh_BASE_(
            self, ax, region_boundary=True, density=6000, usetex=False,
            corlormap='tab10', show_numbering=True, title=True,
            show_boundary_names=True, xlim=None, ylim=None,
            region_boundary_linewidth=1.5, element_linewidth=1,
            element_color='red'
    ):
        """

        :param ax: The axis we are plotting.
        :param region_boundary:
        :param density:
        :param usetex:
        :param corlormap:
        :param show_numbering:
        :param title:
        :param show_boundary_names:
        :param xlim:
        :param ylim:
        :param region_boundary_linewidth:
        :param element_linewidth:
        :param element_color:
        :return:
        """
        if RANK == MASTER_RANK:

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
                    density=4 * 150 * self._mesh_.domain.regions.num, do_plot=False)
            o = np.linspace(0, 1, density)  # plot density
            I_ = np.ones(density)
            # line data ...
            RI = {}  # data for regions internal lines
            for rn in self._mesh_.domain.regions.names:
                RI[rn] = ([], [])  # ([dy_lines], [dx_lines])
                # ____ compute dy lines ___________________________________________
                for x in self._mesh_.elements.spacing[rn][0][:]:
                    RI[rn][0].append(self._mesh_.domain.regions(rn).interpolation.mapping(x * I_, o))
                # ____ compute dx lines ___________________________________________
                for y in self._mesh_.elements.spacing[rn][1][:]:
                    RI[rn][1].append(self._mesh_.domain.regions(rn).interpolation.mapping(o, y * I_))
                # --------------------------------------------------------------------------

            # ... element - internal - mesh ....
            REI = dict()
            nodes = self._f_.space.nodes
            nx, ny = nodes
            if nx[0] == -1:
                nx = nx[1:]
            if nx[-1] == 1:
                nx = nx[:-1]
            if ny[0] == -1:
                ny = ny[1:]
            if ny[-1] == 1:
                ny = ny[:-1]

            # print(nx, ny)
            for rn in self._mesh_.domain.regions.names:
                REI[rn] = ([], [])  # ([dy_lines], [dx_lines])
                for i, x0 in enumerate(self._mesh_.elements.spacing[rn][0][:-1]):
                    x1 = self._mesh_.elements.spacing[rn][0][i+1]
                    for nxi in nx:
                        x = (nxi + 1) * (x1-x0) / 2 + x0
                        REI[rn][0].append(self._mesh_.domain.regions(rn).interpolation.mapping(x * I_, o))
                for i, y0 in enumerate(self._mesh_.elements.spacing[rn][1][:-1]):
                    y1 = self._mesh_.elements.spacing[rn][1][i+1]
                    for nyi in ny:
                        y = (nyi + 1) * (y1-y0) / 2 + y0
                        REI[rn][1].append(self._mesh_.domain.regions(rn).interpolation.mapping(o, y * I_))

            reodb = self._mesh_.domain.regions.edges_on_domain_boundaries
            # text: element numbering data ...
            RBN = None
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
            boundary_name_color_dict: dict = dict()
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
                    boundary_name_color_dict[pb2] = boundary_name_color_dict[pb1]

                    if usetex:
                        pb_text[pb1] = r'\\mathrm{%s}' % pb1 + \
                                       r'\\stackrel{\\mathrm{%s}}{=}' % ptype + '\\mathrm{%s}' % pb2
                        pb_text[pb2] = r'\\mathrm{%s}' % pb2 + \
                                       r'\\stackrel{\\mathrm{%s}}{=}' % ptype + '\\mathrm{%s}' % pb1
                    else:
                        pb_text[pb1] = r'\\mathrm{%s}' % pb1 + \
                                       r'\\genfrac{}{}{0}{}{%s}{=}' % ptype + '\\mathrm{%s}' % pb2
                        pb_text[pb2] = r'\\mathrm{%s}' % pb2 + \
                                       r'\\genfrac{}{}{0}{}{%s}{=}' % ptype + '\\mathrm{%s}' % pb1

            element_center_coordinates = None
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
            # ___________ do the plot ______________________________________________________

            plt.rc('text', usetex=usetex)
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

            for rn in self._mesh_.domain.regions.names:
                for dy_lines in RI[rn][0]:  # plot dy mesh lines
                    ax.plot(dy_lines[0], dy_lines[1], color=element_color, linewidth=element_linewidth)
                for dx_lines in RI[rn][1]:  # plot dx mesh lines
                    ax.plot(dx_lines[0], dx_lines[1], color=element_color, linewidth=element_linewidth)

                for dy_lines in REI[rn][0]:  # plot dy mesh lines
                    ax.plot(dy_lines[0], dy_lines[1], color='gray', linewidth=0.75*element_linewidth)
                for dx_lines in REI[rn][1]:  # plot dx mesh lines
                    ax.plot(dx_lines[0], dx_lines[1], color='gray', linewidth=0.75*element_linewidth)

            for rn in self._mesh_.domain.regions.names:
                for ei in range(4):
                    if reodb[rn][ei] == 1:  # plot the domain boundary
                        bn = self._mesh_.domain.regions.map[rn][ei]
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color=boundary_name_color_dict[bn], linewidth=3)
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1],
                                color='k', linewidth=0.1*3)
                    else:
                        if region_boundary:  # plot the regions boundary
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='b', linewidth=region_boundary_linewidth)

                    if show_boundary_names:
                        if RBN[rn][ei] is None:
                            pass
                        else:
                            bn = self._mesh_.domain.regions.map[rn][ei]
                            if bn in pb_text:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<' + pb_text[bn] + '>$', fontsize=13,
                                        c=boundary_name_color_dict[bn], ha='center', va='center')
                            else:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<$' + bn + '$>$', fontsize=13,
                                        c=boundary_name_color_dict[bn], ha='center', va='center')

            if show_numbering:
                for rn in self._mesh_.domain.regions.names:
                    eccrn = element_center_coordinates[rn]
                    AEGN = self._mesh_.___PRIVATE_generate_ALL_element_global_numbering___()
                    gnrn = AEGN[rn]
                    for i in range(self._mesh_.elements.layout[rn][0]):
                        for j in range(self._mesh_.elements.layout[rn][1]):
                            plt.text(eccrn[0][i, j], eccrn[1][i, j], "@${}$".format(gnrn[i, j]),
                                     color='k', alpha=0.5, fontsize=12, ha='center', va='center')

            if title is True:
                title = f"numbering of {self._f_.orientation}-{self._f_.k}-form: " + \
                        f"'{self._f_.standard_properties.name}'"
                if self._f_.whether.hybrid:
                    title += " (hybrid)"
                plt.title(title)
            elif title is False:
                pass
            else:
                plt.title(title)

            plt.xlabel('$x$')
            plt.ylabel('$y$')
            # plt.tight_layout()
            return cm.get_cmap('Dark2', 8)

    def _matplot__2dCSCG_0Form_Inner_numbering(self, saveto=None, **kwargs):
        """"""
        assert self._f_.k == 0
        nodes = self._f_.space.nodes
        nodes = np.meshgrid(*nodes, indexing='ij')

        elements = self._mesh_.elements
        MAPPING = dict()
        GATHERING = dict()
        for i in elements:
            ele = elements[i]
            MAPPING[i] = ele.coordinate_transformation.mapping(*nodes)
            GATHERING[i] = self._f_.numbering.gathering[i]

        MAPPING = COMM.gather(MAPPING, root=MASTER_RANK)
        GATHERING = COMM.gather(GATHERING, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            mapping = dict()
            gathering = dict()
            for i, MPi in enumerate(MAPPING):
                mapping.update(MPi)
                gathering.update(GATHERING[i])

            local_numbering = self._f_.numbering.local[0]
            assert np.shape(local_numbering) == (self._f_.p[0]+1, self._f_.p[1]+1)
        else:
            return

        fig, ax = plt.subplots(figsize=(15, 9))
        LN_colors = self._matplot_mesh_BASE_(ax, **kwargs)

        # .. now, we attach the numbering of 0-form to the fig.
        is_hybrid = self._f_.whether.hybrid
        for k in mapping:  # go through all elements (kth)
            mpk = mapping[k]
            gtk = gathering[k]
            ck = LN_colors(k % 8)
            for j in range(self._f_.p[1]+1):
                for i in range(self._f_.p[0]+1):
                    local_number = local_numbering[i, j]
                    x = mpk[0][i, j]
                    y = mpk[1][i, j]
                    text = gtk[local_number]
                    if is_hybrid:
                        if i == 0 and j == 0:  # corner UL
                            ax.text(x, y, str(text), c=ck, va='bottom', ha='left')
                        elif i == self._f_.p[0] and j == 0:  # corner DL
                            ax.text(x, y, str(text), c=ck, va='bottom', ha='right')
                        elif i == 0 and j == self._f_.p[1]:  # corner UR
                            ax.text(x, y, str(text), c=ck, va='top', ha='left')
                        elif i == self._f_.p[0] and j == self._f_.p[1]:
                            ax.text(x, y, str(text), c=ck, va='top', ha='right')
                        elif i == 0:
                            ax.text(x, y, str(text), c=ck, ha='left', va='center')
                        elif i == self._f_.p[0]:
                            ax.text(x, y, str(text), c=ck, ha='right', va='center')
                        elif j == 0:
                            ax.text(x, y, str(text), c=ck, va='bottom', ha='center')
                        elif j == self._f_.p[1]:
                            ax.text(x, y, str(text), c=ck, va='top', ha='center')
                        else:
                            ax.text(x, y, str(text), c=ck, va='center', ha='center')
                    else:
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')

        plt.show()
        # ... SAVE TO ...
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
        # ...
        return fig

    def _matplot__2dCSCG_0Form_Outer_numbering(self, **kwargs):
        return self._matplot__2dCSCG_0Form_Inner_numbering(**kwargs)

    def _matplot__2dCSCG_1Form_Inner_numbering(self, saveto=None, **kwargs):
        """"""
        assert self._f_.k == 1 and self._f_.orientation == 'inner'
        nodes = self._f_.space.nodes
        nx, ny = nodes
        cx = (nx[1:] + nx[:-1]) / 2
        dx = np.meshgrid(cx, ny, indexing='ij')
        cy = (ny[1:] + ny[:-1]) / 2
        dy = np.meshgrid(nx, cy, indexing='ij')

        elements = self._mesh_.elements

        MAPPING_dx = dict()
        MAPPING_dy = dict()
        for i in elements:
            ele = elements[i]
            MAPPING_dx[i] = ele.coordinate_transformation.mapping(*dx)
            MAPPING_dy[i] = ele.coordinate_transformation.mapping(*dy)
        MAPPING_dx = COMM.gather(MAPPING_dx, root=MASTER_RANK)
        MAPPING_dy = COMM.gather(MAPPING_dy, root=MASTER_RANK)

        GATHERING = dict()
        for i in elements:
            GATHERING[i] = self._f_.numbering.gathering[i]
        GATHERING = COMM.gather(GATHERING, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            mapping_dx = dict()
            mapping_dy = dict()
            gathering = dict()
            for i, MPxi in enumerate(MAPPING_dx):
                mapping_dx.update(MPxi)
                mapping_dy.update(MAPPING_dy[i])
                gathering.update(GATHERING[i])

            local_numbering_dx, local_numbering_dy = self._f_.numbering.local  # shift for outer-oriented 1-form
            assert np.shape(local_numbering_dx) == (self._f_.p[0], self._f_.p[1] + 1)
            assert np.shape(local_numbering_dy) == (self._f_.p[0] + 1, self._f_.p[1])
        else:
            return

        fig, ax = plt.subplots(figsize=(15, 9))
        LN_colors = self._matplot_mesh_BASE_(ax, **kwargs)

        # .. now, we attach the numbering of inner 1-form to the fig.
        is_hybrid = self._f_.whether.hybrid
        for k in range(self._mesh_.elements.global_num):  # go through all elements (kth)
            mpk_x = mapping_dx[k]
            mpk_y = mapping_dy[k]
            gtk = gathering[k]
            ck = LN_colors(k % 8)
            # ..... dx goes firstly; shift for outer-oriented-1-forms...
            for j in range(self._f_.p[1] + 1):
                for i in range(self._f_.p[0]):
                    local_number = local_numbering_dx[i, j]
                    x = mpk_x[0][i, j]
                    y = mpk_x[1][i, j]
                    text = gtk[local_number]
                    if is_hybrid:
                        if j == 0:
                            ax.text(x, y, str(text), c=ck, va='bottom', ha='center')
                        elif j == self._f_.p[1]:
                            ax.text(x, y, str(text), c=ck, va='top', ha='center')
                        else:
                            ax.text(x, y, str(text), c=ck, va='center', ha='center')
                    else:
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')
            # ..... dy goes secondly; shift for outer-oriented-1-forms...
            for j in range(self._f_.p[1]):
                for i in range(self._f_.p[0]+1):
                    local_number = local_numbering_dy[i, j]
                    x = mpk_y[0][i, j]
                    y = mpk_y[1][i, j]
                    text = gtk[local_number]
                    if is_hybrid:
                        if i == 0:
                            ax.text(x, y, str(text), c=ck, ha='left', va='center')
                        elif i == self._f_.p[0]:
                            ax.text(x, y, str(text), c=ck, ha='right', va='center')
                        else:
                            ax.text(x, y, str(text), c=ck, va='center', ha='center')
                    else:
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')

        plt.show()
        # ... SAVE TO ...
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
        # ...
        return fig

    def _matplot__2dCSCG_1Form_Outer_numbering(self, saveto=None, **kwargs):
        """"""
        assert self._f_.k == 1 and self._f_.orientation == 'outer'
        nodes = self._f_.space.nodes
        nx, ny = nodes
        cx = (nx[1:] + nx[:-1]) / 2
        dx = np.meshgrid(cx, ny, indexing='ij')
        cy = (ny[1:] + ny[:-1]) / 2
        dy = np.meshgrid(nx, cy, indexing='ij')

        elements = self._mesh_.elements

        MAPPING_dx = dict()
        MAPPING_dy = dict()
        for i in elements:
            ele = elements[i]
            MAPPING_dx[i] = ele.coordinate_transformation.mapping(*dx)
            MAPPING_dy[i] = ele.coordinate_transformation.mapping(*dy)
        MAPPING_dx = COMM.gather(MAPPING_dx, root=MASTER_RANK)
        MAPPING_dy = COMM.gather(MAPPING_dy, root=MASTER_RANK)

        GATHERING = dict()
        for i in elements:
            GATHERING[i] = self._f_.numbering.gathering[i]
        GATHERING = COMM.gather(GATHERING, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            mapping_dx = dict()
            mapping_dy = dict()
            gathering = dict()
            for i, MPxi in enumerate(MAPPING_dx):
                mapping_dx.update(MPxi)
                mapping_dy.update(MAPPING_dy[i])
                gathering.update(GATHERING[i])

            local_numbering_dy, local_numbering_dx = self._f_.numbering.local  # shift for inner-oriented 1-form
            assert np.shape(local_numbering_dx) == (self._f_.p[0], self._f_.p[1] + 1)
            assert np.shape(local_numbering_dy) == (self._f_.p[0] + 1, self._f_.p[1])
        else:
            return

        fig, ax = plt.subplots(figsize=(15, 9))
        LN_colors = self._matplot_mesh_BASE_(ax, **kwargs)

        # .. now, we attach the numbering of outer 1-form to the fig.
        is_hybrid = self._f_.whether.hybrid
        for k in range(self._mesh_.elements.global_num):  # go through all elements (kth)
            mpk_x = mapping_dx[k]
            mpk_y = mapping_dy[k]
            gtk = gathering[k]
            ck = LN_colors(k % 8)
            # ..... dy goes firstly; shift for inner-oriented-1-forms...
            for j in range(self._f_.p[1]):
                for i in range(self._f_.p[0]+1):
                    local_number = local_numbering_dy[i, j]
                    x = mpk_y[0][i, j]
                    y = mpk_y[1][i, j]
                    text = gtk[local_number]
                    if is_hybrid:
                        if i == 0:
                            ax.text(x, y, str(text), c=ck, ha='left', va='center')
                        elif i == self._f_.p[0]:
                            ax.text(x, y, str(text), c=ck, ha='right', va='center')
                        else:
                            ax.text(x, y, str(text), c=ck, va='center', ha='center')
                    else:
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')
            # ..... dx goes secondly; shift for inner-oriented-1-forms...
            for j in range(self._f_.p[1]+1):
                for i in range(self._f_.p[0]):
                    local_number = local_numbering_dx[i, j]
                    x = mpk_x[0][i, j]
                    y = mpk_x[1][i, j]
                    text = gtk[local_number]
                    if is_hybrid:
                        if j == 0:
                            ax.text(x, y, str(text), c=ck, va='bottom', ha='center')
                        elif j == self._f_.p[1]:
                            ax.text(x, y, str(text), c=ck, va='top', ha='center')
                        else:
                            ax.text(x, y, str(text), c=ck, va='center', ha='center')
                    else:
                        ax.text(x, y, str(text), c=ck, va='center', ha='center')

        plt.show()
        # ... SAVE TO ...
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
        # ...
        return fig

    def _matplot__2dCSCG_2Form_Inner_numbering(self, saveto=None, **kwargs):
        """"""
        assert self._f_.k == 2
        nodes = self._f_.space.nodes
        nx, ny = nodes
        nx = (nx[1:] + nx[:-1]) / 2
        ny = (ny[1:] + ny[:-1]) / 2
        nodes = np.meshgrid(nx, ny, indexing='ij')

        elements = self._mesh_.elements
        MAPPING = dict()
        GATHERING = dict()
        for i in elements:
            ele = elements[i]
            MAPPING[i] = ele.coordinate_transformation.mapping(*nodes)
            GATHERING[i] = self._f_.numbering.gathering[i]

        MAPPING = COMM.gather(MAPPING, root=MASTER_RANK)
        GATHERING = COMM.gather(GATHERING, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            mapping = dict()
            gathering = dict()
            for i, MPi in enumerate(MAPPING):
                mapping.update(MPi)
                gathering.update(GATHERING[i])

            local_numbering = self._f_.numbering.local[0]
            assert np.shape(local_numbering) == (self._f_.p[0], self._f_.p[1])
        else:
            return

        fig, ax = plt.subplots(figsize=(15, 9))
        LN_colors = self._matplot_mesh_BASE_(ax, **kwargs)

        # .. now, we attach the numbering of 0-form to the fig.
        for k in mapping:  # go through all elements (kth)
            mpk = mapping[k]
            gtk = gathering[k]
            ck = LN_colors(k % 8)
            for j in range(self._f_.p[1]):
                for i in range(self._f_.p[0]):
                    local_number = local_numbering[i, j]
                    x = mpk[0][i, j]
                    y = mpk[1][i, j]
                    text = gtk[local_number]
                    ax.text(x, y, str(text), c=ck, va='center', ha='center')

        plt.show()
        # ... SAVE TO ...
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
        # ...
        return fig

    def _matplot__2dCSCG_2Form_Outer_numbering(self, **kwargs):
        return self._matplot__2dCSCG_2Form_Inner_numbering(**kwargs)

    def local(self, **kwargs):
        """Visualize the local numbering."""
        return self._f_.space.local_numbering.___PRIVATE_matplot___(self._f_.__class__.__name__, **kwargs)
