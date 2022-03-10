


from screws.freeze.main import FrozenOnly


from root.config.main import *
import numpy as np
import matplotlib.pyplot as plt



class _2dCSCG_FormVisualize_Matplot(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, sf):
        self._sf_ = sf
        self._mesh_ = sf.mesh
        self._freeze_self_()

    def __call__(self, **kwargs):
        """"""
        return getattr(self, f'_matplot_{self._sf_.k}form')(**kwargs)


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




    def _matplot_2form(self, **kwargs):
        return self._matplot_0form(**kwargs)

    def _matplot_0form(self, ptype='contourf', **kwargs):
        return getattr(self, '_matplot_0form_' + ptype)(**kwargs)

    def _matplot_0form_contourf(self, density=10000, levels=None, num_levels=20,
        usetex=False, colormap='coolwarm', show_colorbar=True, title=True, show_boundaries=True):
        """"""
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.GLOBAL_num)))
        rs = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xy, v = self._sf_.reconstruct(*rs)

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

            x, y, v = self._mesh_.do.regionwsie_stack(x, y, v)



            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)


            for rn in self._sf_.mesh.domain.regions.names:
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

            if title is True:
                title = self._sf_.orientation + f' {self._sf_.k}form: ' + self._sf_.standard_properties.name
                plt.title(title)

            elif title is False:
                pass
            else:
                plt.title(title)

            plt.xlabel('$x$')
            plt.ylabel('$y$')

            if show_colorbar: plt.colorbar()
            plt.show()
            return fig







    def _matplot_1form(self, ptype='contourf', **kwargs):
        return getattr(self, '_matplot_1form_' + ptype)(**kwargs)

    def _matplot_1form_contourf(self, density=10000, levels_x=None, levels_y=None, num_levels=20,
        usetex=False, colormap='coolwarm', show_colorbar=True, title=True, show_boundaries=True):
        """"""
        density = int(np.ceil(np.sqrt(density / self._mesh_.elements.GLOBAL_num)))
        rs = [np.linspace(-1, 1, density) for _ in range(self._sf_.ndim)]
        xy, v = self._sf_.reconstruct(*rs)

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

            x, y, vx, vy = self._mesh_.do.regionwsie_stack(x, y, vx, vy)

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

            for rn in self._sf_.mesh.domain.regions.names:
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
            if title is True:
                if self._sf_.orientation == 'inner':
                    dfx = r"$(u, \cdot)$ on $\mathrm{d}x$"
                else:
                    dfx = r"$(u, \cdot)$ on $\mathrm{d}y$"
                plt.title(dfx)
            elif title is False:
                pass
            else:
                pass

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

            for rn in self._sf_.mesh.domain.regions.names:
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
            if title is True:
                if self._sf_.orientation == 'inner':
                    dfy = r"$(\cdot, v)$ on $\mathrm{d}y$"
                else:
                    dfy = r"$(\cdot, v)$ on $\mathrm{d}x$"
                plt.title(dfy)
            elif title is False:
                pass
            else:
                pass
            if show_colorbar: plt.colorbar()

            # ----------------------------------------------------------------------------
            if title is True:
                default_title = f'{self._sf_.orientation} {self._sf_.k}-form: ' + \
                                f'{self._sf_.standard_properties.name}'
                plt.suptitle(default_title)
            elif title is False:
                pass
            else:
                pass

            plt.show()
            return fig
