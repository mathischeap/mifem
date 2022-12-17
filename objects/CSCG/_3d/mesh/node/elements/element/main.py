# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.node.elements.element.coordinate_transformation import _3dCSCG_Node_Element_CT
from objects.CSCG._3d.mesh.node.elements.element.whether import _3dCSCG_NodeElement_Whether
from objects.CSCG._3d.mesh.node.elements.element.helpers.parse_boundary_position_indicator import \
    parse_boundary_position_indicator


class _3dCSCG_Node_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        """"""
        self._elements_ = elements
        self._i_ = i
        self._positions_ = elements._locations_[i]
        self._ct_ = None
        self._cp_ = None
        self._ce_ = None
        self._cc_ = None
        self._whether_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def positions(self):
        """This node element is at these positions (all)."""
        return self._positions_

    @property
    def boundary_position_indicator(self):
        """

        Returns
        -------
        indicator :
            If this node-element is an internal-node-element, `indicator` = None.

            Else, we will use function parse_boundary_position_indicator to parse the indicator.

        """
        if not self.whether.on_mesh_boundary:
            return None

        # ----------- this node-element is no mesh-element -------------------------

        indicator = parse_boundary_position_indicator(self)

        return indicator

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Node_Element_CT(self)
        return self._ct_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _3dCSCG_NodeElement_Whether(self)
        return self._whether_

    @property
    def shared_by_mesh_elements(self):
        return self._elements_._shared_by_elements_[self._i_]

    @property
    def on_mesh_boundaries(self):
        return self._elements_._on_mesh_boundaries_[self._i_]

    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this node element."""
        if self._cp_ is None:

            for pos in self.positions:
                element, corner = pos[:-3], pos[-3:]  # before reach boundary names, we must have found it. no worries.
                element = int(element)
                if element in self._elements_._MAP_:
                    self._cp_ = pos
                    self._ce_ = element
                    self._cc_ = corner
                    break
        return self._cp_

    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this node element is a corner of this mesh
        element.
        """
        if self._ce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._ce_

    @property
    def CHARACTERISTIC_corner(self):
        """We mainly consider this node element is such a corner of the
        CHARACTERISTIC_element.
        """
        if self._cc_ is None:
            _ = self.CHARACTERISTIC_position
        return self._cc_

    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this node element is in this region."""
        region = self._elements_._mesh_.do.FIND_region_name_of_element(self.CHARACTERISTIC_element)
        return region
