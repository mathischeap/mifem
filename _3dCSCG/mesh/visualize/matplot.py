

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *

from screws.frozen import FrozenOnly
import numpy as np
import matplotlib.pyplot as plt


class _3dCSCG_Mesh_Visualize_Matplot(FrozenOnly):
    """"""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', " <MeshVisualize> "
        assert mesh.ndim == 3, " <MeshVisualize> "
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.grid(*args, **kwargs)


    def grid(self, elements=None, density=50000, usetex=False,
        saveto = None, linewidth=0.6, aspect='equal',):
        """ We compute the grid from mesh element, so even for periodic boundaries, the grid will be full.

        :param elements: (default: ``None``) If it is ``None``, we plot
            all elements. Otherwise, we plot elements No. ``elements``.
        :param density: How resolved the plot will be.
        :param usetex:
        :param saveto:
        :param linewidth:
        :param aspect:
        :return:
        """
        if elements is None:
            ELE = self._mesh_.elements.indices
        else:
            if not isinstance(elements, (list, tuple)):
                elements = [elements,]

            ELE = list()
            for i in elements:
                if i in self._mesh_.elements.indices:
                    ELE.append(i)

        density = int(np.ceil((
            self._mesh_.ndim*density/self._mesh_.elements.GLOBAL_num
                               )**(1/self._mesh_.ndim)))
        r = np.linspace(-1,1, density)
        O = np.ones(density)
        WB = (r, -1*O, -1*O)
        EB = (r, +1*O, -1*O)
        WF = (r, -1*O, +1*O)
        EF = (r, +1*O, +1*O)
        NB = (-1*O, r, -1*O)
        SB = (+1*O, r, -1*O)
        NF = (-1*O, r, +1*O)
        SF = (+1*O, r, +1*O)
        NW = (-1*O, -1*O, r)
        NE = (-1*O, +1*O, r)
        SW = (+1*O, -1*O, r)
        SE = (+1*O, +1*O, r)

        LINES = dict()
        for i in ELE:
            LINES[i] = dict()
            element = self._mesh_.elements[i]
            LINES[i]['WB'] = element.coordinate_transformation.mapping(*WB)
            LINES[i]['EB'] = element.coordinate_transformation.mapping(*EB)
            LINES[i]['WF'] = element.coordinate_transformation.mapping(*WF)
            LINES[i]['EF'] = element.coordinate_transformation.mapping(*EF)
            LINES[i]['NB'] = element.coordinate_transformation.mapping(*NB)
            LINES[i]['SB'] = element.coordinate_transformation.mapping(*SB)
            LINES[i]['NF'] = element.coordinate_transformation.mapping(*NF)
            LINES[i]['SF'] = element.coordinate_transformation.mapping(*SF)
            LINES[i]['NW'] = element.coordinate_transformation.mapping(*NW)
            LINES[i]['NE'] = element.coordinate_transformation.mapping(*NE)
            LINES[i]['SW'] = element.coordinate_transformation.mapping(*SW)
            LINES[i]['SE'] = element.coordinate_transformation.mapping(*SE)

        LINES = cOmm.gather(LINES, root=mAster_rank)

        if rAnk != mAster_rank: return

        # from now on, we will only in the master core.
        ALL_LINES = dict()
        for lines in LINES:
            ALL_LINES.update(lines)
        del LINES

        plt.rc('text', usetex=usetex)

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

        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for i in ALL_LINES:
            for edge in ALL_LINES[i]:
                plt.plot(*ALL_LINES[i][edge], 'k', linewidth=linewidth)
                if aspect == 'equal':
                    x, y, z = ALL_LINES[i][edge]
                    x_lim.append(np.max(x))
                    x_lim.append(np.min(x))
                    y_lim.append(np.max(y))
                    y_lim.append(np.min(y))
                    z_lim.append(np.max(z))
                    z_lim.append(np.min(z))

                    if len(x_lim) > 20:
                        x_lim = [min(x_lim), max(x_lim)]
                        y_lim = [min(y_lim), max(y_lim)]
                        z_lim = [min(z_lim), max(z_lim)]

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)

        if usetex is False:
            plt.title(self._mesh_.domain.name + ', ID: '+
                      self._mesh_.standard_properties.parameters['ID'] +
                      ', <mesh>')

        #__________ SAVE TO ___________________________________________________________
        plt.tight_layout()
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
        #================================================================================
        return fig

    def element_distribution(self):
        """"""
        return self._mesh_.___PRIVATE_matplot_local_elements___()

    def connection(self, density=5000, usetex=False, saveto=None, aspect='equal'):
        """The connection of node and edge elements.

        :param density:
        :param usetex:
        :param saveto:
        :param aspect:
        :return:
        """
        num_edges = self._mesh_.edge.elements.GLOBAL_num
        density = int(density / num_edges)

        if density <= 1: density = 2

        ep0 = np.linspace(-1,1, density)
        ep1 = np.linspace(-1,1, density)
        ep2 = np.linspace(-1,1, density)

        nodes = self._mesh_.node.elements
        edges = self._mesh_.edge.elements

        nodes_coo = dict() # coordinates
        nodes_omb = nodes._on_mesh_boundaries_
        nodes_opb = nodes._IS_on_periodic_boundary_
        for i in nodes:
            node = nodes[i]
            nodes_coo[i] = node.coordinate_transformation.mapping(from_element='any')

        edges_coo = dict()
        for i in edges:
            edge = edges[i]
            edges_coo[i] = edge.coordinate_transformation.mapping(ep0, ep1, ep2, from_element='any')

        nodes_coo = cOmm.gather(nodes_coo, root=mAster_rank)
        nodes_omb = cOmm.gather(nodes_omb, root=mAster_rank)
        nodes_opb = cOmm.gather(nodes_opb, root=mAster_rank)
        edges_coo = cOmm.gather(edges_coo, root=mAster_rank)

        if rAnk == mAster_rank:
            ___ = dict()
            for _ in nodes_coo:
                ___.update(_)
            nodes_coo = ___

            ___ = dict()
            for _ in nodes_omb:
                ___.update(_)
            nodes_omb = ___

            ___ = dict()
            for _ in nodes_opb:
                ___.update(_)
            nodes_opb = ___

            ___ = dict()
            for _ in edges_coo:
                ___.update(_)
            edges_coo = ___

            assert len(nodes_coo) == nodes.GLOBAL_num, "safety check!"
            assert len(nodes_omb) == nodes.GLOBAL_num, "safety check!"
            assert len(nodes_opb) == nodes.GLOBAL_num, "safety check!"
            assert len(edges_coo) == edges.GLOBAL_num, "safety check!"

        else:
            return

        plt.rc('text', usetex=usetex)
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

        for i in nodes_coo:
            if nodes_omb[i]:
                ax.scatter(*nodes_coo[i], marker='s', color='b')
            elif nodes_opb[i]:
                ax.scatter(*nodes_coo[i], marker='o', color='r')
            else:
                ax.scatter(*nodes_coo[i], marker='s', color='gray')

        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for i in edges_coo:
            ax.plot(*edges_coo[i], color='lightgray')

            if aspect == 'equal':
                x, y, z = edges_coo[i]
                x_lim.append(np.max(x))
                x_lim.append(np.min(x))
                y_lim.append(np.max(y))
                y_lim.append(np.min(y))
                z_lim.append(np.max(z))
                z_lim.append(np.min(z))

                if len(x_lim) > 20:
                    x_lim = [min(x_lim), max(x_lim)]
                    y_lim = [min(y_lim), max(y_lim)]
                    z_lim = [min(z_lim), max(z_lim)]

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)

        if usetex is False:
            plt.title(self._mesh_.domain.name + ', ID: '+
                      self._mesh_.standard_properties.parameters['ID'] +
                      ', <mesh>')

        #__________ SAVE TO ___________________________________________________________
        plt.tight_layout()
        if saveto is not None and saveto != '':
            plt.savefig(saveto, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
        #===============================================================================================
        return fig


if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\mesh\visualize\matplot.py
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
