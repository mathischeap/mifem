# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 3:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK, np, COMM
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm


class mpRfT2_Mesh_rcWds_Scalar_Visualize(FrozenOnly):
    """"""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
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

    def matplot(self, xy,
        plot_type = 'contour',

        title=None,
        levels=None, num_levels=20,
        linewidth=1, linestyle=None,

        usetex=False, colormap='coolwarm',

        show_mesh = False,
        show_colorbar=True,
                colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
                colorbar_labelsize=12.5, colorbar_extend='both',

        ticksize = 12,
        labelsize = 15,

        show_boundaries=True,
        saveto=None, dpi=None,
        ):
        """Plot this scalar on the coordinates `xy`.

        Parameters
        ----------
        xy
        plot_type
        title
        levels
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
        assert xy.__class__.__name__ == 'mpRfT2_Mesh_rcWds_Vector' and \
               xy._mesh_ is self._scalar_._mesh_, \
            f"mesh does not match."

        scalar = self._scalar_
        mesh = scalar._mesh_

        if not scalar._isfull_: raise NotImplementedError()

        xy = xy.rgW
        V = scalar.rgW

        if show_mesh:

            CPD = dict()
            for rp in mesh.rcfc:  # ao through all local cell indices.
                cell = mesh[rp]
                assert cell.___isroot___
                CPD[rp] = cell.coordinate_transformation.___PRIVATE_plot_data___(density=None)
            CPD = COMM.gather(CPD, root=MASTER_RANK)

        if RANK != MASTER_RANK: return

        if show_mesh:
            ___ = dict()
            for _ in CPD:
                ___.update(_)
            CPD = ___

        # -------- levels decider ---------------------------------------------
        if levels is None:
            v = list()
            for rn in xy:
                v.append(V[rn])
            v = np.array(v)
            levels = self.___set_contour_levels___(v, num_levels)
            del v
        else:
            pass
        #----------------------------------------------------------------

        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        if colormap is not None: plt.rcParams['image.cmap'] = colormap
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)

        for rn in mesh.cscg.domain.regions.names:
            x, y = xy[rn]
            if plot_type =='contour':
                plt.contour(x, y, V[rn], levels=levels, linewidths=linewidth, linestyles=linestyle)
            elif plot_type =='contourf':
                plt.contourf(x, y, V[rn], levels=levels)
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

        #-------------------------------------------------------------------------------
        if show_mesh:
            for rp in CPD:
                lines = CPD[rp]
                for xy in lines:
                    plt.plot(*xy, color='k', linewidth=0.4)

        #------------ title ---------------------------------------------------------
        if title is None:
            pass
        else:
            plt.title(title)

        plt.xlabel('$x$', fontsize=labelsize)
        plt.ylabel('$y$', fontsize=labelsize)
        ax.tick_params(labelsize=ticksize)
        #-------------------------------- color bar ---------------------------------
        if show_colorbar:
            mappable = cm.ScalarMappable()
            mappable.set_array(np.array(levels))
            cb = plt.colorbar(mappable, ax=ax,
                              extend=colorbar_extend,
                              aspect=colorbar_aspect,
                              orientation=colorbar_orientation)

            if colorbar_label is not None:
                cb.set_label(colorbar_label, labelpad=10, size=15)

            cb.ax.tick_params(labelsize=colorbar_labelsize)

        #---------------------- save to ---------------------------------------------
        if saveto is None or saveto is False:
            plt.show()
        else:
            if saveto[-4:] == '.df':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')
        plt.close()





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
