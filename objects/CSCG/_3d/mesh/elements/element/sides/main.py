"""

"""

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.elements.element.sides.side.main import _3dCSCG_Mesh_Element_Side



class _3dCSCG_Mesh_Element_Sides(FrozenOnly):
    """The mesh element sides class."""
    def __init__(self, element):
        self._element_ = element
        self._individual_sides_ = dict()
        self._freeze_self_()

    def __iter__(self):
        """only return the indices, as all others."""
        for s in 'NSWEBF':
            yield s

    def __getitem__(self, s):
        assert s in self, f"side name {s} wrong, should be one of 'NSWEBF'."
        if s not in self._individual_sides_:
            self._individual_sides_[s] = _3dCSCG_Mesh_Element_Side(self, s)
        return self._individual_sides_[s]

    def __len__(self):
        return 6

    def __contains__(self, s):
        return s in 'NSWEBF'