# -*- coding: utf-8 -*-



import sys
if './' not in sys.path: sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.edge.elements.element.coordinate_transformation import _3dCSCG_Edge_Element_CT
from objects.CSCG._3d.mesh.edge.elements.element.whether import _3dCSCG_EdgeElement_Whether



class _3dCSCG_Edge_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        """"""
        self._elements_ = elements
        self._i_ = i
        self._positions_ = elements._locations_[i]
        self._direction_ = None
        self._ct_ = None
        self._cp_ = None
        self._ce_ = None
        self._cce_ = None
        self._whether_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def positions(self):
        """This edge element is at these positions."""
        return self._positions_

    @property
    def direction(self):
        """This edge element topologically is along with direction?

        For example, a 'EF' edge is along 'NS' direction.

        :return: one of ('NS', 'WE, 'BF').
        """
        if self._direction_ is None:
            pos = self.positions[0][-2:]
            if pos in ('NB', 'SB', 'NF', 'SF'):
                self._direction_ = self._elements_.___WE___
            elif pos in ('WB', 'EB', 'WF', 'EF'):
                self._direction_ = self._elements_.___NS___
            elif pos in ('NW', 'NE', 'SW', 'SE'):
                self._direction_ = self._elements_.___BF___
            else:
                raise Exception()

        return self._direction_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Edge_Element_CT(self)
        return self._ct_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _3dCSCG_EdgeElement_Whether(self)
        return self._whether_

    @property
    def shared_by_mesh_elements(self):
        return self._elements_._shared_by_elements_[self._i_]

    @property
    def on_mesh_boundaries(self):
        """This edge element is on these mesh boundaries"""
        return self._elements_._on_mesh_boundaries_[self._i_]


    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this edge element."""
        if self._cp_ is None:

            for pos in self.positions:
                element, corner_edge = pos[:-2], pos[-2:]
                # before reach boundary names, we must have found it. So no worries.
                element = int(element)
                if element in self._elements_._MAP_:
                    self._cp_ = pos
                    self._ce_ = element
                    self._cce_ = corner_edge
                    break

        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this edge element is a corner-edge of this mesh
        element."""
        if self._ce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._ce_
    @property
    def CHARACTERISTIC_corner_edge(self):
        """We mainly consider this edge element is such a corner-edge of the
        CHARACTERISTIC_element."""
        if self._cce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._cce_
    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this edge element is in this regions."""
        region = self._elements_._mesh_.do.FIND_region_name_of_element(self.CHARACTERISTIC_element)
        return region







if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\edge\elements\element\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements

    for i in edges:
        edge = edges[i]
        print(edge.i, edge.positions, edge.direction, edge.whether.on_mesh_boundary)