# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/30 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from root.config.main import RANK, SECRETARY_RANK, COMM
from components.freeze.base import FrozenOnly
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm




class _2cCSCG_DV_VisualizeMatplot(FrozenOnly):
    """"""

    def __init__(self, dv):
        """"""
        self._dv_ = dv
        self._mesh_ = dv.mesh
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        return self.quiver(*args, **kwargs)

    def quiver(self, title=None,

        usetex=False,

        axis_on=True, tick_on=True, show_boundaries=False,

        colormap='cool', color=None,

        xlim = None, ylim=None,

        show_colorbar=True,
        colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
        colorbar_labelsize=12.5, colorbar_extend='both',
        colorbar_position = None, colorbar_ticks=None,

        scale=10, scale_units='xy',

        # ticks
        xticks=None, yticks=None,
        tick_size=12, tick_pad=6, minor_tick_length=4, major_tick_length=8,
        label_size = 15,

        saveto=None, dpi=210,

        key_coordinates = (0.5, 0.5),
        key_length = 1,
        key_label = '1'
        ):
        """Could be very badly distributed arrows for non-uniform meshes. Try to use visualization
        of discrete forms.

        Parameters
        ----------
        title
        usetex

        axis_on :
            Show frame.
        tick_on :
            Show ticks.

            Only applies when `axis_on` is True.
        colormap
        color : {None, str}, optional, default: None
            `color` affects when `show_colorbar` is False, and when `color` is not None, we color
            all arrows with it. Otherwise, we color arrows with `colormap`.

        xlim :
            Limit the x-direction.
        ylim :
            Limit the y-direction.
        show_colorbar
        colorbar_label
        colorbar_orientation
        colorbar_aspect
        colorbar_labelsize
        colorbar_extend
        colorbar_position
        colorbar_ticks
        tick_size

        scale : float, optional, default: 10
            Number of data units per arrow length unit, e.g., m/s per plot width; a smaller scale
            parameter makes the arrow longer. Default is None.

            If None, a simple autoscaling algorithm is used, based on the average vector length and the number of vectors.
            The arrow length unit is given by the scale_units parameter.
        scale_units : {'width', 'height', 'dots', 'inches', 'x', 'y', 'xy'}, optional, default: 'xy'
            If the scale kwarg is None, the arrow length unit. Default is None.
            e.g. scale_units is 'inches', scale is 2.0, and (u, v) = (1, 0), then the vector will be 0.5 inches long.

            If scale_units is 'width' or 'height', then the vector will be half the width/height of the axes.

            If scale_units is 'x' then the vector will be 0.5 x-axis units. To plot vectors in the x-y plane,
            with u and v having the same units as x and y, use angles='xy', scale_units='xy', scale=1.

        xticks
        yticks
        tick_size
        tick_pad
        minor_tick_length
        major_tick_length

        label_size
        show_boundaries : bool, optional, default: False

        saveto
        dpi:
            The dpi for pixel based figures.

        key_coordinates:
            Coordinate system and units for X, Y: 'axes' and 'figure' are normalized coordinate
            systems with (0, 0) in the lower left and (1, 1) in the upper right; 'data'
            are the axes data coordinates (used for the locations of the vectors in the
            quiver plot itself); 'inches' is position in the figure in inches, with (0, 0)
            be the lower left corner.
        key_length:
            The length of the key.
        key_label:
            The key label (e.g., length and units of the key).
        Returns
        -------

        """

        xy = self._dv_.coordinates
        v = self._dv_.values
        xy = COMM.gather(xy, root=SECRETARY_RANK)
        v = COMM.gather(v, root=SECRETARY_RANK)

        if RANK != SECRETARY_RANK: return

        X = list()
        Y = list()
        U = list()
        V = list()

        for _xy_, _v_ in zip(xy, v):
            for rn in _xy_:
                _xy = _xy_[rn]
                _v = _v_[rn]
                _x, _y = _xy
                X.append(_x.ravel('F'))
                Y.append(_y.ravel('F'))
                _x, _y = _v
                U.append(_x.ravel('F'))
                V.append(_y.ravel('F'))

        X = np.concatenate(X)
        Y = np.concatenate(Y)
        U = np.concatenate(U)
        V = np.concatenate(V)
        M = np.hypot(U, V)

        #---- check if zero field, if it is quiver will return warning, so we skip it ---------
        U_max, U_min = np.max(U), np.min(U)
        V_max, V_min = np.max(V), np.min(V)
        if U_max - U_min == 0 and  V_max - V_min == 0:
            ZERO_FIELD = True
        else:
            ZERO_FIELD = False

        #----------------------------------------------------------------------------------------
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

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
            mesh = self._mesh_
            RB, RBN, boundary_name_color_dict, pb_text = \
                mesh.visualize.matplot.___PRIVATE_DO_generate_boundary_data___(
                    50, usetex=usetex)[0:4]
            reo_db = mesh.domain.regions.edges_on_domain_boundaries
            for rn in mesh.domain.regions.names:
                for ei in range(4):
                    if reo_db[rn][ei] == 1:
                        # noinspection PyUnresolvedReferences
                        ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k', # not an error, do not use else
                                linewidth=0.75)
        else:
            pass


        #----------------------------------------------------------------------------------------

        if show_colorbar:
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
            norm = matplotlib.colors.Normalize() # M normalized to be [0,1]
            norm.autoscale(M)
            cm = getattr(matplotlib.cm, colormap)
            sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
            sm.set_array([])

            if colorbar_position is not None:
                cb_axes = fig.add_axes(colorbar_position)
                cbar = plt.colorbar(sm, orientation=colorbar_orientation, cax=cb_axes,
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
            pass

        #----------------------------------------------------------------------------------------
        if ZERO_FIELD:
            pass
        else:
            if show_colorbar:
                ax.quiver(X, Y, U, V, color=cm(norm(M)), scale=scale, scale_units=scale_units)

                # Q = ax.quiver(X, Y, U, V, color=cm(norm(M)), scale=scale, scale_units=scale_units)
                # ax.quiverkey(Q, *key_coordinates, key_length, key_label,
                #              color = 'k', # override the color by black
                #              labelpos='E', coordinates='figure')

            else:

                if colormap is not None: plt.rcParams['image.cmap'] = 'Greys'

                if color is None:
                    Q = ax.quiver(X, Y, U, V, M, scale=scale, scale_units=scale_units)
                else:
                    Q = ax.quiver(X, Y, U, V, color='k', scale=scale, scale_units=scale_units)

                ax.quiverkey(Q, *key_coordinates, key_length, key_label,
                             labelpos='E', coordinates='figure')


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
        return plt

    @staticmethod
    def ___set_contour_levels___(v, num_levels, i):
        """ We set the `num_levels` according to the values `v` to be plotted. """
        MAX = list()
        MIN = list()
        for rn in v:
            v_ = v[rn][i]
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
        levels_x=None, levels_y=None, num_levels=20, linewidth=1, linestyles=None,

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

        title_x=True,
        title_y=True,
        saveto=None, dpi=210,
        plot_type = 'contour'):
        """

        Parameters
        ----------
        levels_x
        levels_y
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
        title_x
        title_y
        saveto
        dpi
        plot_type : str, optional
            'contourf' or 'contour'.

        Returns
        -------

        """
        mesh = self._dv_.mesh
        xy = self._dv_.coordinates
        v = self._dv_.values
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
        if levels_x is None:
            levels_x = self.___set_contour_levels___(v, num_levels, i=0)
        else:
            pass
        if levels_y is None:
            levels_y = self.___set_contour_levels___(v, num_levels, i=1)
        else:
            pass

        #------- config plt ----------------------------------------------------------------------1
        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        if colormap is not None: plt.rcParams['image.cmap'] = colormap
        fig = plt.figure(figsize=(12,5.5))

        #---------------------- x-component ------------------------------------------------------1
        ax = plt.subplot(121) #----------- set labels, ticks, frame and so on ----------2
        ax.set_aspect('equal')

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

        #----------------- mesh boundaries --------------------------------------------2
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

        #-------------- contour or contourf plot -----------------------------------2
        for rn in xy:
            if plot_type =='contour':
                if color is None:
                    plt.contour(*xy[rn], v[rn][0], levels=levels_x, linewidths=linewidth, linestyles=linestyles)
                else:
                    plt.contour(*xy[rn], v[rn][0], colors=color, levels=levels_x, linewidths=linewidth, linestyles=linestyles)

            elif plot_type =='contourf':
                VAL = v[rn][0]
                VAL[VAL > levels_x[-1]] = levels_x[-1]
                VAL[VAL < levels_x[0]] = levels_x[0]
                plt.contourf(*xy[rn], VAL, levels=levels_x)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        # --------------- title ----------------------------------------------------2
        if title_x is True:
            title_x =  'x-component-of-' + self._dv_.standard_properties.name
            plt.title(title_x)
        elif title_x is False:
            pass
        else:
            plt.title(title_x)

        #-------------------------------- color bar -------------------------------2
        if (plot_type =='contour' and color is None and show_colorbar) or (plot_type =='contourf' and show_colorbar):
            mappable = cm.ScalarMappable()
            mappable.set_array(np.array(levels_x))
            cb = plt.colorbar(mappable, ax=ax,
                              extend=colorbar_extend,
                              aspect=colorbar_aspect,
                              orientation=colorbar_orientation)

            if colorbar_label is not None:
                cb.set_label(colorbar_label, labelpad=10, size=15)

            cb.ax.tick_params(labelsize=colorbar_labelsize)

        else:
            pass

        #---------------------- y-component ------------------------------------------------------1
        ax = plt.subplot(122) #----------- set labels, ticks, frame and so on ----------2
        ax.set_aspect('equal')

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

        #----------------- mesh boundaries --------------------------------------------2
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

        #-------------- contour or contourf plot -----------------------------------2
        for rn in xy:
            if plot_type =='contour':
                if color is None:
                    plt.contour(*xy[rn], v[rn][1], levels=levels_y, linewidths=linewidth, linestyles=linestyles)
                else:
                    plt.contour(*xy[rn], v[rn][1], colors=color, levels=levels_y, linewidths=linewidth, linestyles=linestyles)

            elif plot_type =='contourf':
                VAL = v[rn][1]
                VAL[VAL > levels_y[-1]] = levels_y[-1]
                VAL[VAL < levels_y[0]] = levels_y[0]
                plt.contourf(*xy[rn], VAL, levels=levels_y)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        # --------------- title ----------------------------------------------------2
        if title_y is True:
            title_y =  'y-component-of-' + self._dv_.standard_properties.name
            plt.title(title_y)
        elif title_y is False:
            pass
        else:
            plt.title(title_y)

        #-------------------------------- color bar -------------------------------2
        if (plot_type =='contour' and color is None and show_colorbar) or (plot_type =='contourf' and show_colorbar):
            mappable = cm.ScalarMappable()
            mappable.set_array(np.array(levels_y))
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


    def contourf(self, *args, **kwargs):
        return self.contour(*args, **kwargs, plot_type='contourf')


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/discrete_fields/vector/visualize/matplot.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('chp2')([[1,2,4,8,8,8,4,2,1],10])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', is_hybrid=True)

    u.TW.func.do.set_func_body_as(ES, 'velocity')
    u.TW.current_time = 0
    u.TW.do.push_all_to_instant()
    u.discretize()

    r = np.linspace(-1,1,5)
    s = np.linspace(-1,1,6)

    dv = u.reconstruct.discrete_vector([r, s])

    mp = dv.visualize.matplot
    # u.visualize.matplot.quiver(show_colorbar=True, colorbar_ticks=[-10, 7,8,9])
    mp.quiver(show_colorbar=True, colorbar_ticks=None, scale=30, key_length=6, key_label=3)