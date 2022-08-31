# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.base import FrozenOnly
from root.config.main import np, sEcretary_rank, cOmm, rAnk
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

# from mpl_toolkits.axes_grid1 import make_axes_locatable

class _2dCSCG_S0F_VIS_Matplot(FrozenOnly):
    """"""
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

    def contourf(self, *args, **kwargs):
        return  self.contour(*args, **kwargs, plot_type='contourf')

    def contour(self, density=10000,

        levels=None, num_levels=20, linewidth=1, linestyle=None,

        usetex=False, colormap='coolwarm',

        show_colorbar=True,
                colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
                colorbar_labelsize=12.5, colorbar_extend='both',

        ticksize = 12,
        labelsize = 15,

        title=True,
        show_boundaries=True,
        saveto=None,
        plot_type = 'contour'):
        """

        :param density:
        :param levels:
        :param num_levels:
        :param linewidth: float
        :param linestyle: {None, 'solid', 'dashed', 'dashdot', 'dotted'}
        :param usetex:
        :param colormap:
        :param show_colorbar: bool
        :param colorbar_label: str
        :param colorbar_orientation: None or {'vertical', 'horizontal'}
            The orientation of the colorbar. It is preferable to set the location of the colorbar,
            as that also determines the orientation; passing incompatible values for location and
            orientation raises an exception.
        :param colorbar_aspect:
        :param colorbar_labelsize:
        :param colorbar_extend: {'neither', 'both', 'min', 'max'}
            If not 'neither', make pointed end(s) for out-of- range values. These are set for a
            given colormap using the colormap set_under and set_over methods.
        :param ticksize:
        :param labelsize:
        :param title:
        :param show_boundaries:
        :param saveto:
        :param plot_type: {'contour', 'contourf'}
        :return:
        """
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.GLOBAL_num)))

        rs = list()
        for _ in range(self._sf_.ndim):
            __ = np.linspace(-1, 1, density+1)
            __ = (__[:-1] + __[1:]) / 2
            rs.append(__)

        xy, v = self._sf_.reconstruct(*rs)

        xy = cOmm.gather(xy, root=sEcretary_rank)
        v = cOmm.gather(v, root=sEcretary_rank)

        if rAnk != sEcretary_rank:
            return

        XY = dict()
        VV = dict()
        for _ in xy: XY.update(_)
        for _ in v: VV.update(_)

        x = list()
        y = list()
        v = list()
        for i in XY:
            x.append(XY[i][0])
            y.append(XY[i][1])
            v.append(VV[i][0])
        v = np.array(v)
        if levels is None:
            levels = self.___set_contour_levels___(v, num_levels)
        else:
            pass

        x, y, v = self._mesh_.do.regionwsie_stack(x, y, v)

        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        if colormap is not None: plt.rcParams['image.cmap'] = colormap
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        #------- label and  ticks -------------------------------------------------------
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        plt.xlabel('$x$', fontsize=labelsize)
        plt.ylabel('$y$', fontsize=labelsize)
        ax.tick_params(labelsize=ticksize)

        #-------------- plot -------------------------------------------------------------
        for rn in self._sf_.mesh.domain.regions.names:
            if plot_type =='contour':
                plt.contour(x[rn], y[rn], v[rn], levels=levels, linewidths=linewidth, linestyles=linestyle)
            elif plot_type =='contourf':
                VAL = v[rn]
                VAL[VAL > levels[-1]] = levels[-1]
                VAL[VAL < levels[0]] = levels[0]
                plt.contourf(x[rn], y[rn], VAL, levels=levels)
            else:
                raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

        #-------- boundaries ------------------------------------------------------------------
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
                            linewidth=0.75)

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
        # --------------- title -------------------------------------------------------------------
        if title is True:
            title = self._sf_.orientation + f' {self._sf_.k}form: ' + self._sf_.standard_properties.name
            plt.title(title)

        elif title is False:
            pass
        else:
            plt.title(title)
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
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')
        plt.close()
        #--------------------------------------------------------------------------

        return fig



if __name__ == '__main__':
    # mpiexec -n 3 python _2dCSCG\forms\standard\_0_form\base\visualize\matplot.py

    from objects.CSCG._2d.master import MeshGenerator, ExactSolutionSelector, SpaceInvoker, FormCaller
    from numpy import pi
    mesh = MeshGenerator('crazy_periodic', bounds=[[0, 2*pi], [0, 2*pi]], c=0.1)([15, 15])
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)], show_info=True)
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)("Euler:shear_layer_rollup", show_info=True)
    w = FC('0-f-o', is_hybrid=False, name='vorticity')
    w.TW.func.do.set_func_body_as(es, 'vorticity')
    w.TW.current_time = 0
    w.TW.do.push_all_to_instant()
    w.discretize()
    w.visualize.matplot.contourf(levels =[-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6], usetex=False)


    mesh = MeshGenerator('rectangle', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=False)
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-i', is_hybrid=True)
    f0.TW.func.do.set_func_body_as(ES, 'potential')
    f0.TW.current_time = 0
    f0.TW.do.push_all_to_instant()
    f0.discretize()
    f0.visualize.matplot.contourf(usetex=True)