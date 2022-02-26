
import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
# import tecplot as tp
# from tecplot.constant import PlotType
from screws.frozen import FrozenOnly

from _3dCSCG.mesh.visualize.matplot import _3dCSCG_Mesh_Visualize_Matplot


class _3dCSCG_Mesh_Visualize(FrozenOnly):
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', " <MeshVisualize> "
        assert mesh.ndim == 3, " <MeshVisualize> "
        self._mesh_ = mesh
        self._matplot_ = _3dCSCG_Mesh_Visualize_Matplot(mesh)
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    def tecplot(self, elements=None, density=100000):
        """
        Tecplot the mesh; the elements. Therefore, we do not plot the internal grids
        for basis functions which are not a part of the mesh.

        In Tecplot: layout -> page -> frame -> dataset -> zones. And each frame can
        also have multiple plots, each time, only one plot will be activated by doing
        `plot.activate()`.

        Parameters
        ----------
        elements : optional
            The elements to be plotted. If it is None (default), we quick plot all
            elements with element boundaries being planes.

            if `elements='all'`, we carefully (slowly) plot all elements: all element
            boundaries will be correctly plotted. Therefore, curvilinear elements will
            appear as curvilinear elements.
        density :

        """
        if rAnk != mAster_rank: return
        if elements is None:
            self._mesh_.domain.visualize.tecplot(nodes=self._mesh_._element_spacing_)

        else: # clearly, below code does not work!

            raise NotImplementedError(density)

            # #______ get data for all elements _________________________________________
            # density = int(np.ceil((self._mesh_.ndim*density/self._mesh_.elements.num
            #                        )**(1/self._mesh_.ndim)))
            # rst = [np.linspace(-1, 1, density) for _ in range(self._mesh_.ndim)]
            # r, s, t = np.meshgrid(*rst, indexing='ij')
            # self._mesh_.evaluate_coordinate_transformation_at(r, s, t)
            # xyz = self._mesh_.coordinate_transformation.mapping
            # self._mesh_.coordinate_transformation._reset_()
            # #_______ which elements to be plotted ______________________________________
            # elements = [i for i in range(self._mesh_.elements.num)] \
            #     if elements == 'all' else elements
            # #_______ lets tecplot _____________________________________________________
            # tp.session.connect()
            # tp.new_layout()
            # page = tp.active_page()
            # page.name = 'Page_Elements'
            # frame = page.active_frame()
            # frame.name = 'Frame_Elements'
            # dataset = frame.create_dataset('Elements')
            # variable_names = 'xyz'
            # for i in range(self._mesh_.ndim):
            #     dataset.add_variable(variable_names[i])
            # for n in elements:
            #     zone = dataset.add_ordered_zone('Element#'+str(n), (density, density, density))
            #     for i in range(self._mesh_.ndim):
            #         zone.values(variable_names[i])[:] = xyz[i][n].ravel('F')
            # frame.plot_type = getattr(PlotType, 'Cartesian'+str(self._mesh_.ndim)+'D')
            # plot = frame.plot()
            # plot.show_shade = False
            # plot.show_edge = True
            # plot.show_mesh = False
            # for j in range(len(elements)):
            #     surfaces = plot.fieldmap(j).surfaces
            #     surfaces.surfaces_to_plot = True
            # plot.use_translucency  = False
            # plot.use_lighting_effect = False
            # plot.view.fit()
            # frame.plot().activate()

    @property
    def matplot(self):
        return self._matplot_







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\mesh\visualize\main.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 3, 3]
    mesh = MeshGenerator('crazy', c=0.25)(elements)
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0, 3], [0, 3], [0, 3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    # mesh = MeshGenerator('psc')(elements)
    mesh.visualize.matplot.grid()
    mesh.visualize.matplot.connection()
    mesh.visualize.matplot.element_distribution()
    mesh.domain.visualize()
    mesh.domain.regions.visualize()