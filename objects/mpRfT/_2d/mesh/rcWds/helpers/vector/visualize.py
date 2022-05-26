# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 7:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, mAster_rank, np, cOmm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm


class mpRfT2_Mesh_rcWds_Vector_Visualize(FrozenOnly):
    """"""

    def __init__(self, vector):
        """"""
        self._vector_ = vector
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)


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



    def matplot(self, xy, plot_type='contourf',  **kwargs):
        """"""
        assert xy.__class__.__name__ == 'mpRfT2_Mesh_rcWds_Vector' and \
               xy._mesh_ is self._vector_._mesh_, \
            f"mesh does not match."

        if plot_type == 'quiver':
            return self.___Pr_quiver___(xy, **kwargs)
        elif plot_type in ('contourf', 'contour'):
            return self.___Pr_contour_f___(xy, plot_type=plot_type, **kwargs)
        else:
            raise NotImplementedError(f"not implemented for plot_type={plot_type}.")



    def ___Pr_contour_f___(self, xy,
        plot_type='contourf',
        title_x=None, title_y=None, suptitle=None,
        levels_x=None, levels_y=None,
        num_levels=20,
        linewidth=1, linestyle=None,

        usetex=False, colormap='coolwarm',

        show_mesh = False,
        show_colorbar=True,
        colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
        colorbar_labelsize=12.5, colorbar_extend='both',

        ticksize=12,
        labelsize=15,

        show_boundaries=True,
        saveto=None, dpi=None,
        ):
        """

        Parameters
        ----------
        xy
        plot_type
        title_x
        title_y
        suptitle
        levels_x
        levels_y
        num_levels
        linewidth
        linestyle
        usetex
        colormap
        show_mesh
        show_colorbar
        colorbar_label
        colorbar_orientation
        colorbar_aspect
        colorbar_labelsize
        colorbar_extend
        ticksize
        labelsize
        show_boundaries
        saveto
        dpi

        Returns
        -------

        """
        vector = self._vector_
        mesh = vector._mesh_

        if not vector._isfull_: raise NotImplementedError()

        xy = xy.rgW
        V = vector.rgW


        if show_mesh:

            CPD = dict()
            for rp in mesh.rcfc:  # ao through all local cell indices.
                cell = mesh[rp]
                assert cell.___isroot___
                CPD[rp] = cell.coordinate_transformation.___PRIVATE_plot_data___(density=None)
            CPD = cOmm.gather(CPD, root=mAster_rank)

        if rAnk != mAster_rank: return

        if show_mesh:
            ___ = dict()
            for _ in CPD:
                ___.update(_)
            CPD = ___

        # -------- levels decider ---------------------------------------------
        if levels_x is None:
            v = list()
            for rn in xy:
                v.append(V[rn][0])
            v = np.array(v)
            levels_x = self.___set_contour_levels___(v, num_levels)
            del v
        else:
            pass
        if levels_y is None:
            v = list()
            for rn in xy:
                v.append(V[rn][1])
            v = np.array(v)
            levels_y = self.___set_contour_levels___(v, num_levels)
            del v
        else:
            pass
        #===============================================================================
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        if colormap is not None: plt.rcParams['image.cmap'] = colormap
        fig = plt.figure(figsize=(12,5.5))

        #------------- x  component ----------------------------------------------------------------
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

        for rn in mesh.cscg.domain.regions.names:
            x, y = xy[rn]
            if plot_type =='contour':
                plt.contour(x, y, V[rn][0], levels=levels_x, linewidths=linewidth, linestyles=linestyle)
            elif plot_type =='contourf':
                plt.contourf(x, y, V[rn][0], levels=levels_x)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        RB, RBN, boundary_name_color_dict, pb_text = \
            mesh.cscg.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]
        reo_db = mesh.cscg.domain.regions.edges_on_domain_boundaries

        for rn in mesh.cscg.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = mesh.cscg.domain.regions.map[rn][ei]
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
                        bn = mesh.cscg.domain.regions.map[rn][ei]
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
        #-----------------------------------------------------------------------------------------
        if show_mesh:
            for ind in CPD:
                lines = CPD[ind]
                for __ in lines:
                    plt.plot(*__, color='k', linewidth=0.4)

        #------------ title ---------------------------------------------------------
        if title_x is None:
            pass
        else:
            plt.title(title_x)

        plt.xlabel('$x$', fontsize=labelsize)
        plt.ylabel('$y$', fontsize=labelsize)
        ax.tick_params(labelsize=ticksize)
        #-------------------------------- color bar ---------------------------------
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

        #------------- y  component ---------------------------------------------------------------
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

        for rn in mesh.cscg.domain.regions.names:
            x, y = xy[rn]
            if plot_type =='contour':
                plt.contour(x, y, V[rn][1], levels=levels_y, linewidths=linewidth, linestyles=linestyle)
            elif plot_type =='contourf':
                plt.contourf(x, y, V[rn][1], levels=levels_y)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        RB, RBN, boundary_name_color_dict, pb_text = \
            mesh.cscg.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]
        reo_db = mesh.cscg.domain.regions.edges_on_domain_boundaries

        for rn in mesh.cscg.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = mesh.cscg.domain.regions.map[rn][ei]
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
                        bn = mesh.cscg.domain.regions.map[rn][ei]
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
        #-----------------------------------------------------------------------------------------
        if show_mesh:
            for ind in CPD:
                lines = CPD[ind]
                for __ in lines:
                    plt.plot(*__, color='k', linewidth=0.4)

        #------------ title ---------------------------------------------------------------
        if title_y is None:
            pass
        else:
            plt.title(title_y)

        plt.xlabel('$x$', fontsize=labelsize)
        plt.ylabel('$y$', fontsize=labelsize)
        ax.tick_params(labelsize=ticksize)
        #-------------------------------- color bar -----------------------------------------
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

        #=========== super title ==============================================================
        if suptitle is None:
            pass
        else:
            plt.suptitle(suptitle)

        #---------------------- save to --------------------------------------------------------
        if saveto is None or saveto is False:
            plt.show()
        else:
            if saveto[-4:] == 'pdf':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')

        plt.close(fig)

    def ___Pr_quiver___(self, xy, title=None,

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
        """

        Parameters
        ----------
        xy
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

        Returns
        -------

        """
        vector = self._vector_
        mesh = vector._mesh_

        if not vector._isfull_: raise NotImplementedError()

        xy = xy.rgW
        v = vector.rgW

        if rAnk != mAster_rank: return

        #____ preparing data __________________________________________________________
        U, V, X, Y = [], [], [], []
        for rn in mesh.cscg.domain.regions.names:
            U.append(v[rn][0])
            V.append(v[rn][1])
            X.append(xy[rn][0])
            Y.append(xy[rn][1])

        U = np.array(U).ravel()
        V = np.array(V).ravel()
        X = np.array(X).ravel()
        Y = np.array(Y).ravel()
        M = np.hypot(U, V)

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
            mesh.cscg.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                50, usetex=usetex)[0:4]
        reo_db = mesh.cscg.domain.regions.edges_on_domain_boundaries

        for rn in mesh.cscg.domain.regions.names:
            for ei in range(4):
                if reo_db[rn][ei] == 1:
                    bn = mesh.cscg.domain.regions.map[rn][ei]
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
                        bn = mesh.cscg.domain.regions.map[rn][ei]
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

        if show_colorbar:
            # M = M / np.max(M) normalize to max == 1
            norm = matplotlib.colors.Normalize()
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
        if title is None:
            pass
        else:
            plt.title(title)

        #---------------------- save to --------------------------------------------------------
        if saveto is None or saveto is False:
            plt.show()
        else:
            if saveto[-4:] == 'pdf':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')

        plt.close(fig)









if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
