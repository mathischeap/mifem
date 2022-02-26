"""


"""



import sys
if './' not in sys.path: sys.path.append('../')

from screws.frozen import FrozenOnly





class _3dCSCG_Mesh_Boundary(FrozenOnly):
    def __init__(self, bdrs, name):
        self._bdrs_ = bdrs
        self._name_ = name
        self._freeze_self_()

    @property
    def element_sides(self):
        """This mesh boundary covers these local mesh element sides."""
        return self._bdrs_.RANGE_element_sides[self._name_]

    @property
    def trace_elements(self):
        """This mesh boundary covers these local trace elements."""
        return self._bdrs_.RANGE_trace_elements[self._name_]




