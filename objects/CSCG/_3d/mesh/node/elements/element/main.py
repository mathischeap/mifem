


from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.node.elements.element.coordinate_transformation import _3dCSCG_Node_Element_CT




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
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def positions(self):
        """This node element is at these positions."""
        return self._positions_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Node_Element_CT(self)
        return self._ct_

    @property
    def shared_by_mesh_elements(self):
        return self._elements_._shared_by_elements_[self._i_]

    @property
    def on_mesh_boundaries(self):
        return self._elements_._on_mesh_boundaries_[self._i_]

    @property
    def IS_on_mesh_boundary(self):
        return True if len(self.on_mesh_boundaries) > 0 else False

    @property
    def IS_on_periodic_boundary(self):
        return self._elements_._IS_on_periodic_boundary_[self._i_]


    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this node element."""
        if self._cp_ is None:

            for pos in self.positions:
                element, corner = pos[:-3], pos[-3:] # before reach boundary names, we must have found it. So no worries.
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
        element."""
        if self._ce_ is None:
            _ = self.CHARACTERISTIC_position
        return self._ce_
    @property
    def CHARACTERISTIC_corner(self):
        """We mainly consider this node element is such a corner of the
        CHARACTERISTIC_element."""
        if self._cc_ is None:
            _ = self.CHARACTERISTIC_position
        return self._cc_
    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this node element is in this regions."""
        region = self._elements_._mesh_.do.FIND_region_name_of_element(self.CHARACTERISTIC_element)
        return region
