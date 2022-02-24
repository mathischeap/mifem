# -*- coding: utf-8 -*-
"""
CONTINUOUS FORMS are special, they do not inherit the main class of form!

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from abc import ABC

from root.config import *
import matplotlib.pyplot as plt
from SCREWS.frozen import FrozenOnly

from INHERITING.CSCG.field.main_BASE import CSCG_Continuous_FORM_BASE



class _2dCSCG_Continuous_FORM_BASE(CSCG_Continuous_FORM_BASE, ABC):
    """"""
    def __init_subclass__(cls, ndim=2):
        super().__init_subclass__(ndim=ndim)

    def __init__(self, mesh, ftype, valid_time):
        """"""
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', "Need a 2dCSCG mesh."
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_field')
        super().__init__(mesh, ftype, valid_time)

        self._visualize_ = _2dCSCG_Field_Visualize(self)







class _2dCSCG_Field_Visualize(FrozenOnly):
    def __init__(self, cf):
        self._cf_ = cf
        self._default_ = 'matplot'
        self._mesh_ = self._cf_.mesh
        self._freeze_self_()

    def __call__(self, **kwargs):
        return getattr(self, self._default_)(**kwargs)

    def matplot(self, **kwargs):
        """ """
        # print(self.___2dCSCG_ScalarField)
        return getattr(self, '_plot_' + self._cf_.__class__.__name__)(**kwargs)



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



    def _plot__2dCSCG_ScalarField(self, ptype='contourf', **kwargs):
        return getattr(self, '_plot__2dCSCG_ScalarField' + f'_{ptype}')(**kwargs)

    def _plot__2dCSCG_ScalarField_contourf(self, density=10000, usetex=False, colormap='coolwarm',
        show_colorbar=True, levels=None, num_levels=20, title=True,
        show_boundaries=True, domain_boundary_linewidth=3, boundary_name_fontsize=12,
        minor_tick_length=0, major_tick_length=0, tick_pad=5, tick_size=12):
        """

        :param density:
        :param usetex:
        :param colormap:
        :param show_colorbar:
        :param levels:
        :param num_levels:
        :param title:
        :param show_boundaries:
        :param domain_boundary_linewidth:
        :param boundary_name_fontsize:
        :param minor_tick_length:
        :param major_tick_length:
        :param tick_pad:
        :param tick_size:
        :return:
        """
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.GLOBAL_num)))
        rs = [np.linspace(-1, 1, density) for _ in range(self._cf_.ndim)]
        xy, v = self._cf_.reconstruct(*rs)

        xy = cOmm.gather(xy, root=sEcretary_rank)
        v = cOmm.gather(v, root=sEcretary_rank)

        if rAnk != sEcretary_rank:
            pass
        else:
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

            x, y, v = self._mesh_.DO.regionwsie_stack(x, y, v)

            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)

            for rn in self._cf_.mesh.domain.regions.names:
                plt.contourf(x[rn], y[rn], v[rn], levels=levels)

            if show_boundaries:
                RB, RBN, boundary_name_color_dict, pb_text = \
                    self._mesh_.visualize.___PRIVATE_DO_generate_boundary_data___(
                        50, usetex=usetex)[0:4]

                reo_db = self._mesh_.domain.regions.edges_on_domain_boundaries

                for rn in self._mesh_.domain.regions.names:
                    for ei in range(4):
                        if reo_db[rn][ei] == 1:
                            bn = self._mesh_.domain.regions.map[rn][ei]
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                    linewidth=domain_boundary_linewidth)
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                    linewidth=0.25*domain_boundary_linewidth)

                        if RBN[rn][ei] is None:
                            pass
                        else:
                            bn = self._mesh_.domain.regions.map[rn][ei]
                            if bn in pb_text:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<' + pb_text[bn] + '>$', fontsize=boundary_name_fontsize,
                                        c=boundary_name_color_dict[bn], ha='center', va='center')
                            else:
                                # noinspection PyUnresolvedReferences
                                ax.text(RBN[rn][ei][0], RBN[rn][ei][1],
                                        '$<$' + bn + '$>$', fontsize=boundary_name_fontsize,
                                        c=boundary_name_color_dict[bn], ha='center', va='center')

            plt.tick_params(which='both', labeltop=False, labelright=False, top=False, right=False)
            if minor_tick_length != 0:
                plt.tick_params(axis='both', which='minor', direction='in', length=minor_tick_length)
            if major_tick_length != 0:
                plt.tick_params(axis='both', which='major', direction='in', length=major_tick_length)
            plt.tick_params(axis='both', which='both', labelsize=tick_size)
            plt.tick_params(axis='x', which='both', pad=tick_pad)
            plt.tick_params(axis='y', which='both', pad=tick_pad)

            plt.xlabel('$x$')
            plt.ylabel('$y$')
            if title is True: plt.title("continuous scalar field")
            if show_colorbar: plt.colorbar()
            plt.show()
            return fig


    def _plot__2dCSCG_VectorField(self, ptype='contourf', **kwargs):
        return getattr(self, '_plot__2dCSCG_VectorField' + f'_{ptype}')(**kwargs)

    def _plot__2dCSCG_VectorField_contourf(self, density=10000, usetex=False, colormap='coolwarm',
        show_colorbar=True, levels_x=None, levels_y=None, num_levels=20, show_boundaries=True, title=True):
        """

        :param density:
        :param usetex:
        :param colormap:
        :param show_colorbar:
        :param levels_x:
        :param levels_y:
        :param num_levels:
        :param show_boundaries:
        :param title:
        :return:
        """
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.GLOBAL_num)))
        rs = [np.linspace(-1, 1, density) for _ in range(self._cf_.ndim)]
        xy, v = self._cf_.reconstruct(*rs)

        xy = cOmm.gather(xy, root=sEcretary_rank)
        v = cOmm.gather(v, root=sEcretary_rank)

        if rAnk != sEcretary_rank:
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

            x, y, vx, vy = self._mesh_.DO.regionwsie_stack(x, y, vx, vy)

            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap
            fig = plt.figure(figsize=(12,5.5))

            # ---------------- x component ---------------------------------------------------
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

            for rn in self._cf_.mesh.domain.regions.names:
                plt.contourf(x[rn], y[rn], vx[rn], levels=levels_x)

            if show_boundaries:
                RB, RBN, boundary_name_color_dict, pb_text = \
                    self._mesh_.visualize.___PRIVATE_DO_generate_boundary_data___(
                        50, usetex=usetex)[0:4]

                reo_db = self._mesh_.domain.regions.edges_on_domain_boundaries

                for rn in self._mesh_.domain.regions.names:
                    for ei in range(4):
                        if reo_db[rn][ei] == 1:
                            bn = self._mesh_.domain.regions.map[rn][ei]
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                    linewidth=3)
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                    linewidth=0.25*3)

                        if RBN[rn][ei] is None:
                            pass
                        else:
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
            if title is True: plt.title(r"$(u, \cdot)$")
            if show_colorbar: plt.colorbar()

            # ---------------- y component ---------------------------------------
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

            for rn in self._cf_.mesh.domain.regions.names:
                plt.contourf(x[rn], y[rn], vy[rn], levels=levels_y)

            if show_boundaries:
                for rn in self._mesh_.domain.regions.names:
                    for ei in range(4):
                        if reo_db[rn][ei] == 1:
                            bn = self._mesh_.domain.regions.map[rn][ei]
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color=boundary_name_color_dict[bn],
                                    linewidth=3)
                            # noinspection PyUnresolvedReferences
                            ax.plot(RB[rn][ei][0], RB[rn][ei][1], color='k',
                                    linewidth=0.25*3)

                        if RBN[rn][ei] is None:
                            pass
                        else:
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
            if title is True: plt.title(r"$(\cdot, v)$")
            if show_colorbar: plt.colorbar()

            # ----------------------------------------------------------------------------
            if title is True: plt.suptitle('continuous vector field')
            plt.show()
            return fig