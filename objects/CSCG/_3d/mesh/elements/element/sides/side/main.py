


from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.elements.element.sides.side.coordinate_transformation import _3dCSCG_Mesh_Element_Side_CT



class _3dCSCG_Mesh_Element_Side(FrozenOnly):
    """The mesh element side class."""
    def __init__(self, sides, position):
        self._sides_ = sides
        self._element_ = sides._element_
        self._mesh_ = self._element_._mesh_
        self._position_ = str(self._element_.i) + position
        self._side_index_= 'NSWEBF'.index(position)
        self._ct_ = None
        self._freeze_self_()

    @property
    def position(self):
        """This side is on the `position` of the element."""
        return self._position_

    @property
    def side_index(self):
        """The index of position in 'NSWEBF'."""
        return self._side_index_

    @property
    def side_name(self):
        """The index of position in 'NSWEBF'."""
        return self.position[-1]

    @property
    def trace_element(self):
        """Return the number of the trace element this element side
        represents.
        """
        return self._mesh_.trace.elements.map[
            self._element_.i][self.side_index]

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Mesh_Element_Side_CT(self)
        return self._ct_