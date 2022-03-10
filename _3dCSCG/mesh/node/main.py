# -*- coding: utf-8 -*-
"""
The edge elements of a mesh.
"""


import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly

from _3dCSCG.mesh.node.elements.main import _3dCSCG_Node_Elements



class _3dCSCG_Node(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = None
        self._freeze_self_()


    @property
    def elements(self):
        """The edge elements. Only generate them when called first time."""
        if self._elements_ is None:
            self._elements_ = _3dCSCG_Node_Elements(self)
        return self._elements_








if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\node.py
    from _3dCSCG.main import MeshGenerator
    elements = [2, 2, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    nodes = mesh.node.elements

    # print(rAnk, mesh.___local_periodic_element_sides___)

    print(nodes.GLOBAL_num)

    # for i in nodes:
    #     node = nodes[i]
        # print(i, edge.i, edge.positions, edge.CHARACTERISTIC_element, edge.CHARACTERISTIC_position)
        # print(i, node.i, node.positions, node.CHARACTERISTIC_position, node.CHARACTERISTIC_element, node.CHARACTERISTIC_corner, node.CHARACTERISTIC_region)
