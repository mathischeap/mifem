# -*- coding: utf-8 -*-
"""The edge elements of a mesh."""


import sys
if './' not in sys.path:
    sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.edge.elements.main import _3dCSCG_Edge_Elements


class _3dCSCG_Edge(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = None
        self._freeze_self_()

    @property
    def elements(self):
        """The edge elements. Only generate them when called first time."""
        if self._elements_ is None:
            self._elements_ = _3dCSCG_Edge_Elements(self)
        return self._elements_


if __name__ == '__main__':
    # mpiexec -n 12 python objects\CSCG\_3d\mesh\edge\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edge = mesh.edge
