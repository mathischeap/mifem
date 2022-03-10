"""

"""


from screws.freeze.main import FrozenOnly




class _3dCSCG_Mesh_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._elements_._mesh_.do.find.slave_of_element(i)