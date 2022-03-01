



from screws.frozen import FrozenOnly




class _3dCSCG_Mesh_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._elements_.___DO_find_slave_of_element___(i)