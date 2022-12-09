# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from root.config.main import np, SECRETARY_RANK, COMM, RANK
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

class _2dCSCG_S1F_VIS_Matplot(FrozenOnly):
    """Mesh-element-wise plotter. May be not good for non-uniform meshes."""
    def __init__(self, sf):
        self._sf_ = sf
        self._mesh_ = sf.mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.contourf(*args, **kwargs)

    @staticmethod
    def ___set_contour_levels___(v, num_levels):
        """ We set the `num_levels` according to the values `v` to be plotted. """
        MINv = np.min(v)
        MAXv = np.max(v)
        if MINv == MAXv:
            if MINv == 0:
                MAXv = 0.1
            else:
                if MINv > 0:
                    MAXv = MINv * 1.01
                else:
                    MAXv = MINv * 0.99
            num_levels = 2
        levels = np.linspace(MINv, MAXv, num_levels)
        return levels

    def contourf(self, density=10000, num_levels=20,
        usetex=False, colormap='coolwarm',

        levels_x=None, levels_y=None,
        show_colorbar=True,
                 colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
                 colorbar_labelsize=12.5, colorbar_extend='both',

        title_x=True,
        title_y=True,
        suptitle = False,
        show_boundaries=True,
        saveto = None,
        plot_type='contourf',
        ):
        """

        Parameters
        ----------
        density
        levels_x
        levels_y
        num_levels
        usetex
        colormap
        show_colorbar
        colorbar_label
        colorbar_orientation
        colorbar_aspect
        colorbar_labelsize
        colorbar_extend
        title_x
        title_y
        suptitle
        show_boundaries
        saveto
        plot_type

        Returns
        -------

        """
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.global_num)))

        rs = list()
        for _ in range(self._sf_.ndim):
            __ = np.linspace(-1, 1, density+1)
            __ = (__[:-1] + __[1:]) / 2
            rs.append(__)

        xy, v = self._sf_.reconstruct(*rs)

        xy = COMM.gather(xy, root=SECRETARY_RANK)
        v = COMM.gather(v, root=SECRETARY_RANK)

        if RANK != SECRETARY_RANK:
            pass
        else:
            XY = dict()
            VV = dict()
            for _ in xy: XY.update(_)
            for _ in v: VV.update(_)

            x = list()
            y = list()
            vx = list()
            vy = list()
            for i in XY:
                x.append(XY[i][0])
                y.append(XY[i][1])
                vx.append(VV[i][0])
                vy.append(VV[i][1])
            vx = np.array(vx)
            vy = np.array(vy)
            if levels_x is None:
                levels_x = self.___set_contour_levels___(vx, num_levels)
            else:
                pass
            if levels_y is None:
                levels_y = self.___set_contour_levels___(vy, num_levels)
            else:
                pass

            x, y, vx, vy = self._mesh_.do.regionwsie_stack(x, y, vx, vy)

            if saveto is not None: matplotlib.use('Agg')
            plt.rc('text', usetex=usetex)
            plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
            if colormap is not None: plt.rcParams['image.cmap'] = colormap
            fig = plt.figure(figsize=(12,4))

            plotter = getattr(plt, plot_type)
            # ---------------- x component ---------------------------------------------------------
            ax = plt.subplot(121)
            plt.axis("equal")
            # noinspection PyUnresolvedReferences
            ax.spines['top'].set_visible(False)
            # noinspection PyUnresolvedReferences
            ax.spines['right'].set_visible(False)
            # noinspection PyUnresolvedReferences
            ax.spines['left'].set_visible(True)
            # noinspection PyUnresolvedReferences
            ax.spines['bottom'].set_visible(True)

            for rn in self._sf_.mesh.domain.regions.names:
                plotter(x[rn], y[rn], vx[rn], levels=levels_x)


            RB, RBN, boundary_name_color_dict, pb_text = \
                self._mesh_.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                    50, usetex=usetex)[0:4]

            reo_db = self._mesh_.domain.regions.edges_on_domain_boundaries

            for rn in self._mesh_.domain.regions.names:
                for ei in range(4):
                    if reo_db[rn][ei] == 1:
                        bn = self._mesh_.domain.regions.map[rn][ei]
                        if show_boundaries:
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                    linewidth=3)
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                linewidth=0.25*3)

                    if RBN[rn][ei] is None:
                        pass
                    else:
                        if show_boundaries:
                            bn = self._mesh_.domain.regions.map[rn][ei]
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
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            if title_x is True:
                if self._sf_.orientation == 'inner':
                    dfx = r"$(u, \cdot)$ on $\mathrm{d}x$"
                else:
                    dfx = r"$(u, \cdot)$ on $\mathrm{d}y$"
                plt.title(dfx)
            elif title_x is False:
                pass
            else:
                plt.title(title_x)

            if show_colorbar:
                mappable = cm.ScalarMappable()
                mappable.set_array(np.array(levels_x))
                cb = plt.colorbar(mappable, ax=ax,
                                  extend=colorbar_extend,
                                  aspect=colorbar_aspect,
                                  orientation=colorbar_orientation)

                if colorbar_label is not None:
                    cb.set_label(colorbar_label, labelpad=10, size=15)

                cb.ax.tick_params(labelsize=colorbar_labelsize)

            # ---------------- y component ---------------------------------------------------------
            ax = plt.subplot(122)
            plt.axis("equal")
            # noinspection PyUnresolvedReferences
            ax.spines['top'].set_visible(False)
            # noinspection PyUnresolvedReferences
            ax.spines['right'].set_visible(False)
            # noinspection PyUnresolvedReferences
            ax.spines['left'].set_visible(True)
            # noinspection PyUnresolvedReferences
            ax.spines['bottom'].set_visible(True)

            for rn in self._sf_.mesh.domain.regions.names:
                plotter(x[rn], y[rn], vy[rn], levels=levels_y)

            for rn in self._mesh_.domain.regions.names:
                for ei in range(4):
                    if reo_db[rn][ei] == 1:
                        bn = self._mesh_.domain.regions.map[rn][ei]
                        # noinspection PyUnresolvedReferences
                        if show_boundaries:
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                    linewidth=3)
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                linewidth=0.25*3)

                    if RBN[rn][ei] is None:
                        pass
                    else:
                        if show_boundaries:
                            bn = self._mesh_.domain.regions.map[rn][ei]
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
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            if title_y is True:
                if self._sf_.orientation == 'inner':
                    dfy = r"$(\cdot, v)$ on $\mathrm{d}y$"
                else:
                    dfy = r"$(\cdot, v)$ on $\mathrm{d}x$"
                plt.title(dfy)
            elif title_y is False:
                pass
            else:
                plt.title(title_y)

            if show_colorbar:
                mappable = cm.ScalarMappable()
                mappable.set_array(np.array(levels_y))
                cb = plt.colorbar(mappable, ax=ax,
                                  extend=colorbar_extend,
                                  aspect=colorbar_aspect,
                                  orientation=colorbar_orientation)

                if colorbar_label is not None:
                    cb.set_label(colorbar_label, labelpad=10, size=15)

                cb.ax.tick_params(labelsize=colorbar_labelsize)

            # --------------------------------------------------------------------------------------
            if suptitle is True:
                default_title = f'{self._sf_.orientation} {self._sf_.k}-form: ' + \
                                f'{self._sf_.standard_properties.name}'
                plt.suptitle(default_title)
            elif suptitle is False:
                pass
            else:
                plt.suptitle(suptitle)
            #---------------------- save to --------------------------------------------------------
            if saveto is None or saveto == '':
                plt.show()
            else:
                plt.savefig(saveto, bbox_inches='tight')

            plt.close()

            return fig

    def contour(self, **kwargs):
        return self.contourf(**kwargs, plot_type='contour')

    def quiver(self, density=100, title=None,

        usetex=False, colormap='cool', xlim = None, ylim=None,

        show_colorbar=True,
        colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
        colorbar_labelsize=12.5, colorbar_extend='both',
        colorbar_position = None, colorbar_ticks=None,

        quiverkey='1<->1',

        ticksize=12,
        labelsize=15,
        show_boundaries=True,
        saveto=None, dpi=None,
        ):
        """Could be very badly distributed arrows for non-uniform meshes. Try to use visualization
        of discrete forms.

        Parameters
        ----------
        density
        title
        usetex
        colormap
        xlim
        ylim
        show_colorbar
        colorbar_label
        colorbar_orientation
        colorbar_aspect
        colorbar_labelsize
        colorbar_extend
        colorbar_position
        colorbar_ticks
        ticksize

        quiverkey : str
            This defines the indicator: an arrow.
            For example:
                `quiverkey = '1 <-> text'`

                This gives a showcase arrow whose length is 1. Then we can compare it
                with the arrows to see what are they length.

                The text 'text' will be added beside the showcase arrow.

        labelsize
        show_boundaries
        saveto
        dpi:
            The dpi for pixel based figures.

        Returns
        -------

        """
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.global_num)))

        mesh = self._mesh_

        rs = list()
        for _ in range(self._sf_.ndim):
            __ = np.linspace(-1, 1, density+1)
            __ = (__[:-1] + __[1:]) / 2
            rs.append(__)

        xy, v = self._sf_.reconstruct(*rs)

        xy = COMM.gather(xy, root=SECRETARY_RANK)
        v = COMM.gather(v, root=SECRETARY_RANK)

        if RANK != SECRETARY_RANK: return

        XY = dict()
        VV = dict()
        for _ in xy: XY.update(_)
        for _ in v: VV.update(_)

        X = list()
        Y = list()
        U = list()
        V = list()
        for i in XY:
            X.append(XY[i][0])
            Y.append(XY[i][1])
            U.append(VV[i][0])
            V.append(VV[i][1])

        U = np.array(U).ravel()
        V = np.array(V).ravel()
        X = np.array(X).ravel()
        Y = np.array(Y).ravel()
        M = np.hypot(U, V)


        #---- check if zero field, if it is quiver will return warning, so we skip it ---------
        U_max, U_min = np.max(U), np.min(U)
        V_max, V_min = np.max(V), np.min(V)
        if U_max - U_min == 0 and  V_max - V_min == 0:
            ZERO_FIELD = True
        else:
            ZERO_FIELD = False

        #-------------------------------------------------------------------------
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        # noinspection PyUnresolvedReferences
        ax.spines['top'].set_visible(False)
        # noinspection PyUnresolvedReferences
        ax.spines['right'].set_visible(False)
        # noinspection PyUnresolvedReferences
        ax.spines['left'].set_visible(True)
        # noinspection PyUnresolvedReferences
        ax.spines['bottom'].set_visible(True)
        ax.set_xlabel(r"$x$", fontsize=labelsize)
        ax.set_ylabel(r"$y$", fontsize=labelsize)
        plt.tick_params(axis='both', which='both', labelsize=ticksize)
        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)

        RB, RBN, boundary_name_color_dict, pb_text = \
            mesh.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]
        reo_db = mesh.domain.regions.edges_on_domain_boundaries

        for rn in mesh.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = mesh.domain.regions.map[rn][ei]
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
                        bn = mesh.domain.regions.map[rn][ei]
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

        if ZERO_FIELD:
            pass
        else:
            if show_colorbar:
                # M = M / np.max(M) normalize to max == 1
                norm = matplotlib.colors.Normalize()
                if  colorbar_ticks is None:
                    pass
                else:
                    tMin, tMax = min(colorbar_ticks), max(colorbar_ticks)
                    assert tMax > tMin, f"colorbar_ticks={colorbar_ticks} wrong!"
                    assert tMin >=0 , f"quiver tick can not be lower than 0!"
                    LARGE = M > tMax
                    M[LARGE] = tMax
                    LOW = M < tMin
                    M[LOW] = tMin
                    M = np.concatenate((M, [tMin, tMax]))

                norm.autoscale(M)
                cm = getattr(matplotlib.cm, colormap)
                sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
                sm.set_array([])
                ax.quiver(X, Y, U, V, color=cm(norm(M)))

                if colorbar_position is not None:
                    cbaxes = fig.add_axes(colorbar_position)
                    cbar = plt.colorbar(sm, orientation=colorbar_orientation, cax=cbaxes,
                                  extend=colorbar_extend,
                                  aspect=colorbar_aspect,)
                else:
                    cbar = plt.colorbar(sm, orientation=colorbar_orientation,
                                  extend=colorbar_extend,
                                  aspect=colorbar_aspect,)

                if colorbar_label is not None:
                    colorbar_label.set_label(colorbar_label, labelpad=10, size=15)

                if colorbar_ticks is not None: cbar.set_ticks(colorbar_ticks)
                cbar.ax.tick_params(labelsize=colorbar_labelsize)

            else:
                if colormap is not None: plt.rcParams['image.cmap'] = 'Greys'
                Q = ax.quiver(X, Y, U, V, M, color='k')

                assert '<->' in quiverkey, " <Quiver> : quiverkey={} format wrong.".format(quiverkey)
                value, quivertext = quiverkey.split('<->')
                try:
                    value = int(value)
                except ValueError:
                    raise Exception(
                            " <Quiver>: quiverkey={} format wrong. "
                            "value (before '<->') is not int.".format(quiverkey))

                ax.quiverkey(Q, 0.8, 0.9, value, quivertext, labelpos='E', coordinates='figure')

        # =========== super title ==============================================================
        if title is None or title == '':
            pass
        else:
            plt.title(title)

        #---------------------- save to --------------------------------------------------------
        if saveto is None or saveto is False:
            plt.show()
        else:
            if saveto[-4:] == '.pdf':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')

        plt.close(fig)