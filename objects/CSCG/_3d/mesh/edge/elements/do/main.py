# -*- coding: utf-8 -*-
"""
"""

import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from root.config.main import COMM, np, MASTER_RANK, RANK, MPI
import matplotlib.pyplot as plt

from objects.CSCG._3d.mesh.edge.elements.do.find.main import _3dCSCG_Edge_Elements_DO_FIND




class _3dCSCG_Edge_Elements_DO(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._FIND_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._FIND_ is None:
            self._FIND_ = _3dCSCG_Edge_Elements_DO_FIND(self._elements_)
        return self._FIND_

    def illustrate_element(self, i):
        """

        Parameters
        ----------
        i

        Returns
        -------

        """
        zoom = 0.9
        if i in self._elements_:

            LS = np.linspace(-1,1,20)*zoom

            edge = self._elements_[i]
            positions = edge.positions
            DATA = dict()
            edge_LOCATION = list()
            TraceElement_Data = dict()

            on_mesh_boundary = edge.IS.on_mesh_boundary

            TE_coordinates = (
                [-1, 0, 0],
                [ 1, 0, 0],
                [0, -1, 0],
                [0,  1, 0],
                [0, 0, -1],
                [0, 0,  1],
            )

            for pos in positions:
                if pos[0] in '0123456789': # at the corner of a mesh element
                    element_num = int(pos[:-2])
                    if element_num in self._elements_._mesh_.elements:
                        element = self._elements_._mesh_.elements[element_num]
                        data = element.do.generate_element_plot_data(zoom=zoom)
                        DATA[element_num] = (data, pos[-2:])

                        edge_LOCATION .append(
                            edge.coordinate_transformation.mapping(
                            LS, LS, LS,
                            from_element=element_num,
                            corner_edge=pos[-2:])
                        )

                        TE_MAP = self._elements_._mesh_.trace.elements.map[element_num]

                        TraceElement_Data[element_num] = list()

                        IND = ['NSWEBF'.index(_) for _ in pos[-2:]]

                        for _, t in enumerate(TE_MAP):
                            if _ in IND:
                                coo = element.coordinate_transformation.mapping(*TE_coordinates[_])
                                TraceElement_Data[element_num].append([t, coo])

        else:
            DATA = dict()
            edge_LOCATION = None
            on_mesh_boundary = False
            TraceElement_Data = dict()


        DATA = COMM.gather(DATA, root=MASTER_RANK)
        TraceElement_Data = COMM.gather(TraceElement_Data, root=MASTER_RANK)
        edge_LOCATION = COMM.gather(edge_LOCATION, root=MASTER_RANK)
        on_mesh_boundary = COMM.reduce(on_mesh_boundary, root=MASTER_RANK, op=MPI.LOR)

        nL = list()
        if RANK == MASTER_RANK:
            ___ = dict()
            for data in DATA:
                ___.update(data)
            DATA = ___

            ___ = dict()
            for data in TraceElement_Data:
                ___.update(data)
            TraceElement_Data = ___

            for __ in edge_LOCATION:
                if __ is not None:
                    for _ in __:
                        nL.append(_)
        else:
            return

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

        for nl in nL:
            ax.plot(*nl, color='b', linewidth=1.25)

        for e in DATA:
            data, pos = DATA[e]
            X = data[:,0,:]
            Y = data[:,1,:]
            Z = data[:,2,:]
            x_lim.append(np.max(X))
            x_lim.append(np.min(X))
            y_lim.append(np.max(Y))
            y_lim.append(np.min(Y))
            z_lim.append(np.max(Z))
            z_lim.append(np.min(Z))

            for line in data:
                ax.plot(*line, color='gray', linewidth=0.75)

            ax.text(np.mean(X), np.mean(Y), np.mean(Z), str(e),
                        bbox=dict(boxstyle="square",
                                   ec=(0., 1, 1),
                                   fc=(0., 1, 1),
                                   ),
                    ha='center', va='center', ma='center')

            Ts_COOs = TraceElement_Data[e]
            for t_coo in Ts_COOs:
                t, coo = t_coo
                cx, cy, cz = coo
                cx = cx[0]
                cy = cy[0]
                cz = cz[0]
                ax.text(cx, cy, cz, str(t), color='g',
                        bbox=dict(boxstyle="square",
                                   ec=(1., 0.5, 0.5, 0.5),
                                   fc=(1., 0.8, 0.8, 0.25),
                                   ),
                        ha='center', va='center', ma='center')

        ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.set_xlabel(r'$x$', fontsize=12)
        ax.set_ylabel(r'$y$', fontsize=12)
        ax.set_zlabel(r'$z$', fontsize=12)
        if on_mesh_boundary:
            plt.title(f'Edge-element #{i} on mesh boundary.')
        else:
            plt.title(f'Edge-element #{i}')
        plt.show()
        plt.close(fig)








if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\edge\elements\do\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements

    edges.do.illustrate_element(0)

    # for i in range(edges.GLOBAL_num):
    #     edges.do.illustrate_element(i)