# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/30/2022 11:50 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from root.config.main import RANK, SECRETARY_RANK, COMM
from screws.freeze.base import FrozenOnly
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm




class _2cCSCG_DS_VisualizeMatplot(FrozenOnly):
    """"""

    def __init__(self, ds):
        """"""
        self._ds_ = ds
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.contour(*args, **kwargs)

    @staticmethod
    def ___set_contour_levels___(v, num_levels):
        """ We set the `num_levels` according to the values `v` to be plotted. """
        MAX = list()
        MIN = list()
        for rn in v:
            v_ = v[rn][0]
            MAX.append(np.max(v_))
            MIN.append(np.min(v_))

        MINv = np.min(MIN)
        MAXv = np.max(MAX)
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

    def contour(self,
        levels=None, num_levels=20, linewidth=1, linestyles=None,

        usetex=False,

        colormap='coolwarm', color=None, # color will override color map, only works in contour.

        axis_on=True, tick_on=True, show_boundaries=False,

        xlim=None, ylim=None,

        show_colorbar=True,
                colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
                colorbar_labelsize=12.5, colorbar_extend='both',

        # ticks
        xticks = None, yticks = None,
        tick_size=12, tick_pad=6, minor_tick_length=4, major_tick_length=8,

        label_size = 15,

        title=True,
        saveto=None, dpi=210,
        plot_type = 'contour'):
        """

        Parameters
        ----------
        levels
        num_levels :
            Only applies when `levels` is None.

        linewidth : {float, int, array-like}
            Only applies to contour.

            The line width of the contour lines.

            If a number, all levels will be plotted with this linewidth.

            If a sequence, the levels in ascending order will be plotted with the linewidths in
            the order specified.

            If None, this falls back to rcParams["lines.linewidth"] (default: 1).

        linestyles
            Only applies to contour.

            If linestyles is None, the default is 'solid' unless the lines are monochrome.
            In that case, negative contours will take their linestyle from
            rcParams["contour.negative_linestyle"] (default: 'dashed') setting.

            linestyles can also be an iterable of the above strings specifying a
            set of linestyles to be used. If this iterable is shorter than the number of
            contour levels it will be repeated as necessary.

        usetex : bool, optional, default: False
            If we use LaTeX render?

        colormap
        color : (None, str), optional
            Only applies to contour.

            When  it is not None, it overrides the colormap.

        axis_on :
            Show frame.
        tick_on :
            Show ticks.

            Only applies when `axis_on` is True.

        show_boundaries
        xlim
        ylim
        show_colorbar
        colorbar_label
        colorbar_orientation
        colorbar_aspect
        colorbar_labelsize
        colorbar_extend
        xticks
        yticks
        tick_size
        tick_pad
        minor_tick_length
        major_tick_length
        label_size
        title
        saveto
        dpi
        plot_type : str, optional
            'contourf' or 'contour'.

        Returns
        -------

        """
        mesh = self._ds_.mesh
        xy = self._ds_.coordinates
        v = self._ds_.values
        xy = COMM.gather(xy, root=SECRETARY_RANK)
        v = COMM.gather(v, root=SECRETARY_RANK)

        if RANK != SECRETARY_RANK: return
        # ----------- prepare  data --------------------------------------------------------------1
        _ = dict()
        for __ in xy:
            _.update(__)
        xy = _

        _ = dict()
        for __ in v:
            _.update(__)
        v = _

        #------------ prepare levels -------------------------------------------------------------1
        if levels is None:
            levels = self.___set_contour_levels___(v, num_levels)
        else:
            pass

        #------- config plt ----------------------------------------------------------------------1
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        if colormap is not None: plt.rcParams['image.cmap'] = colormap
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        #----------- set labels, ticks, frame and so on ------------------------------------------1
        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)
        if axis_on:
            if tick_on:
                if xticks is not None: plt.xticks(xticks)
                if yticks is not None: plt.yticks(yticks)
                plt.tick_params(which='both', labeltop=False, labelright=False, top=False, right=False)
                plt.tick_params(axis='both', which='minor', direction='out', length=minor_tick_length)
                plt.tick_params(axis='both', which='major', direction='out', length=major_tick_length)
                plt.tick_params(axis='both', which='both', labelsize=tick_size)
                plt.tick_params(axis='x', which='both', pad=tick_pad)
                plt.tick_params(axis='y', which='both', pad=tick_pad)
                plt.xlabel('$x$', fontsize=label_size)
                plt.ylabel('$y$', fontsize=label_size)
            else:
                plt.tick_params(which='both',
                    labeltop=False, labelright=False, labelleft=False, labelbottom=False,
                    top=False, right=False, left=False, bottom=False)
        else:
            ax.set_axis_off()

        #----------------- mesh boundaries -------------------------------------------------------1
        if show_boundaries:
            RB, RBN, boundary_name_color_dict, pb_text = \
                mesh.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                    50, usetex=usetex)[0:4]
            reo_db = mesh.domain.regions.edges_on_domain_boundaries
            for rn in mesh.domain.regions.names:
                for ei in range(4):
                    if reo_db[rn][ei] == 1:
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                linewidth=0.75)
        else:
            pass

        #-------------- contour or contourf plot -------------------------------------------------1
        for rn in mesh.domain.regions.names:
            if plot_type =='contour':
                if color is None:
                    plt.contour(*xy[rn], v[rn][0], levels=levels, linewidths=linewidth, linestyles=linestyles)
                else:
                    plt.contour(*xy[rn], v[rn][0], colors=color, levels=levels, linewidths=linewidth, linestyles=linestyles)

            elif plot_type =='contourf':
                VAL = v[rn][0]
                VAL[VAL > levels[-1]] = levels[-1]
                VAL[VAL < levels[0]] = levels[0]
                plt.contourf(*xy[rn], VAL, levels=levels)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        # --------------- title ------------------------------------------------------------------1
        if title is True:
            title =  self._ds_.standard_properties.name
            plt.title(title)
        elif title is False:
            pass
        else:
            plt.title(title)

        #-------------------------------- color bar ----------------------------------------------1
        if (plot_type =='contour' and color is None and show_colorbar) or (plot_type =='contourf' and show_colorbar):
            mappable = cm.ScalarMappable()
            mappable.set_array(np.array(levels))
            cb = plt.colorbar(mappable, ax=ax,
                              extend=colorbar_extend,
                              aspect=colorbar_aspect,
                              orientation=colorbar_orientation)

            if colorbar_label is not None:
                cb.set_label(colorbar_label, labelpad=10, size=15)

            cb.ax.tick_params(labelsize=colorbar_labelsize)

        else:
            pass

        #---------------------- save to ----------------------------------------------------------1
        if saveto is None:
            plt.show()
        else:
            if saveto[-4:] == '.pdf':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.savefig(saveto, dpi=dpi, bbox_inches='tight')
        plt.close()

        #=========================================================================================1
        return fig


    def contourf(self, *are, **kwargs):
        """"""
        return self.contour(*are, **kwargs, plot_type='contourf')


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/discrete_fields/scalar/visualize/matplot.py
    pass