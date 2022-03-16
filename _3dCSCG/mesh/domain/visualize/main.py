"""

"""
from root.config.main import *
import tecplot as tp
from tecplot.constant import PlotType
from screws.freeze.main import FrozenOnly
import matplotlib.pyplot as plt
from matplotlib import cm


class _3dCSCG_Domain_Visualize(FrozenOnly):
    def __init__(self, domain):
        """ """
        assert domain.__class__.__name__ == '_3dCSCG_Domain', " <DomainPlot> "
        assert domain.ndim == 3, " <DomainPlot> "
        self._domain_ = domain
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    def tecplot(self, nodes=None):
        """
        Here we call the tecplot API and plot self in Tecplot 360. In Tecplot:
        layout -> page -> frame -> dataset -> zones. And each frame can also
        have multiple plots, each time, only one plot will be activated by
        doing plot.activate().

        Parameters
        ----------
        nodes : int or None, optional
            A positive int to determine how good we are going to follow the real shape
            of each region. When it is None (default), it will be set to be the optimal
            according the regions type.

        """
        if rAnk != mAster_rank: return
        tp.session.connect()
        tp.new_layout()
        page = tp.active_page()
        page.name = 'Domain: ' + self._domain_.name
        frame = page.active_frame()
        dataset = frame.create_dataset('Domain')
        dataset.add_variable('x')
        dataset.add_variable('y')
        dataset.add_variable('z')

        for region_name in self._domain_.regions():
            region = self._domain_.regions(region_name)
            if nodes is None: nodes = 5

            if isinstance(nodes, int):
                assert nodes >= 2, " <Region3D> : density={} wrong".format(nodes)
                r = s = t = np.linspace(0, 1, nodes)
            else:
                r, s, t = nodes[region_name]

            size_x = np.size(r)
            size_y = np.size(s)
            size_z = np.size(t)
            r, s, t = np.meshgrid(r, s, t, indexing='ij')
            r = r.ravel('F')
            s = s.ravel('F')
            t = t.ravel('F')
            x, y, z = region.interpolation(r, s, t)
            zone = dataset.add_ordered_zone('Zone: ' + region.name,
                                            (size_x, size_y, size_z))
            zone.values('x')[:] = x
            zone.values('y')[:] = y
            zone.values('z')[:] = z

        frame.plot_type = PlotType.Cartesian3D
        plot = frame.plot()
        plot.show_shade = True
        plot.show_edge = True
        plot.use_translucency = True
        for i in range(self._domain_._num_regions_):
            surfaces = plot.fieldmap(i).surfaces
            surfaces.surfaces_to_plot = True
        plot.show_mesh = False
        plot.view.fit()
        frame.plot().activate()


    def matplot(self, density=1e4, corlormap='tab10',
        show_internal_region_sides=False,
        show_boundary_names=True,
        distinguish_boundaries_by_color=True,
        aspect='equal',):
        """"""
        # we can do everything in the master core.
        if rAnk != mAster_rank: return

        density = int(np.ceil(np.sqrt(density/(self._domain_.regions.num*6))))
        rrr = sss = np.linspace(0,1,density)
        rrr, sss = np.meshgrid(rrr, sss, indexing='ij')

        regions = self._domain_.regions
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

        boundaries_num = self._domain_.boundaries.num
        boundaries_name = self._domain_.boundaries.names
        if boundaries_num > 10 and corlormap=='tab10': corlormap = 'viridis'
        color = cm.get_cmap(corlormap, boundaries_num)
        colors = []
        boundary_name_color_dict = dict()
        boundary_name_color_solid_dict = dict()
        for j in range(boundaries_num):
            temp_color = list(color(j))
            temp_color[-1] = 0.45
            colors.append(temp_color)
        for j, bn in enumerate(boundaries_name):
            boundary_name_color_dict[bn] = colors[j]
            boundary_name_color_solid_dict[bn] = color(j)

        # initializing for aspect
        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for rn in self._domain_.regions:
            region = self._domain_.regions[rn]
            for i, sn in enumerate('NSWEBF'):
                is_boundary = regions.sides_on_domain_boundaries[rn][i]

                rs = region.sides[sn]

                xyz = rs.coordinate_transformation.mapping(rrr, sss)

                # find the ration.
                if aspect == 'equal' or show_boundary_names:
                    x, y, z = xyz
                    x_lim.append(np.max(x))
                    x_lim.append(np.min(x))
                    y_lim.append(np.max(y))
                    y_lim.append(np.min(y))
                    z_lim.append(np.max(z))
                    z_lim.append(np.min(z))
                else:
                    x, y, z = None, None, None

                if is_boundary:
                    if distinguish_boundaries_by_color:
                        bn = regions.map[rn][i]
                        ax.plot_surface(*xyz, color=boundary_name_color_dict[bn])
                    else:
                        ax.plot_surface(*xyz, color=(1,1,1,0.6))

                    if show_boundary_names:
                        x_range = [np.min(x), np.max(x)]
                        y_range = [np.min(y), np.max(y)]
                        z_range = [np.min(z), np.max(z)]
                        x_mid, y_mid, z_mid = np.mean(x_range), np.mean(y_range), np.mean(z_range)

                        bn = regions.map[rn][i]
                        ax.text(x_mid, y_mid, z_mid, bn,
                                color=boundary_name_color_solid_dict[bn],
                                ha='center', va='center', ma='center')

                else:
                    if show_internal_region_sides:
                        ax.plot_surface(*xyz, color=(1,1,1,0.3))

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        plt.title(self._domain_.name + ', ID: '+
                  self._domain_.parameters['ID'] +
                  ', <domain>')
        fig.tight_layout()

        plt.show()
        plt.close('all')