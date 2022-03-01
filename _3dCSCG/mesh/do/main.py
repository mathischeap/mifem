

from screws.frozen import FrozenOnly
from _3dCSCG.mesh.do.find import _3dCSCG_Mesh_DO_FIND


class _3dCSCG_Mesh_DO(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._FIND_ = _3dCSCG_Mesh_DO_FIND(self)
        self._freeze_self_()

    def reset_cache(self):
        self._mesh_.___PRIVATE_reset_cache___()

    def parse_element_side_pair(self, eP):
        return self._mesh_.___DO_parse_element_side_pair___(eP)

    @property
    def find(self):
        return self._FIND_


    def regionwsie_stack(self, *args):
        return self._mesh_.___DO_regionwsie_stack___(*args)