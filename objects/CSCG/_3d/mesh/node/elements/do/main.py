

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm, np, MPI
import matplotlib.pyplot as plt

from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.mesh.node.elements.do.find import _3dCSCG_NodeMesh_DoFind

class _3dCSCG_NodeMesh_Do(FrozenOnly):
    """"""
    def __init__(self, elements):
        self._elements_ = elements
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3dCSCG_NodeMesh_DoFind(self._elements_)
        return self._find_

    def illustrate_element(self, i):
        """We use this method to illustrate the node element #i.

        This method has to be called at the elements level because it may use mesh-elements
        at different cores.

        Parameters
        ----------
        i : int

        Returns
        -------

        """

        if i in self._elements_:
            node = self._elements_[i]
            positions = node.positions
            DATA = dict()
            node_LOCATION = list()
            TraceElement_Data = dict()

            on_mesh_boundary = node.IS.on_mesh_boundary

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
                    element_num = int(pos[:-3])
                    if element_num in self._elements_._mesh_.elements:
                        element = self._elements_._mesh_.elements[element_num]
                        data = element.do.generate_element_edge_data(zoom=0.9)
                        DATA[element_num] = (data, pos[-3:])

                        node_LOCATION .append(
                            node.coordinate_transformation.mapping(
                            from_element=element_num,
                                corner=pos[-3:])
                        )

                        TE_MAP = self._elements_._mesh_.trace.elements.map[element_num]

                        TraceElement_Data[element_num] = list()

                        IND = ['NSWEBF'.index(_) for _ in pos[-3:]]

                        for _, t in enumerate(TE_MAP):
                            if _ in IND:
                                coo = element.coordinate_transformation.mapping(*TE_coordinates[_])
                                TraceElement_Data[element_num].append([t, coo])

        else:
            DATA = dict()
            node_LOCATION = None
            on_mesh_boundary = None
            TraceElement_Data = dict()

        DATA = cOmm.gather(DATA, root=mAster_rank)
        TraceElement_Data = cOmm.gather(TraceElement_Data, root=mAster_rank)
        node_LOCATION = cOmm.gather(node_LOCATION, root=mAster_rank)
        on_mesh_boundary = cOmm.reduce(on_mesh_boundary, root=mAster_rank, op=MPI.LOR)

        nL = list()
        if rAnk == mAster_rank:
            ___ = dict()
            for data in DATA:
                ___.update(data)
            DATA = ___

            ___ = dict()
            for data in TraceElement_Data:
                ___.update(data)
            TraceElement_Data = ___

            for __ in node_LOCATION:
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
            ax.scatter(*nl, marker='s', color='b')

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
            plt.title(f'Node-element #{i} on mesh boundary.')
        else:
            plt.title(f'Node-element #{i}')
        plt.show()
        plt.close(fig)






if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\node\elements\do\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 3, 3]
    mesh = MeshGenerator('crazy', c=0., bounds=([0,3], [0,3], [0,3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    nodes = mesh.node.elements

    for i in range(mesh.node.elements.GLOBAL_num):
        nodes.do.illustrate_element(i)

