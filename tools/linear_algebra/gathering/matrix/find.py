
from screws.freeze.main import FrozenOnly
from functools import lru_cache





class Gathering_Matrix_FIND(FrozenOnly):
    """"""
    def __init__(self, GM):
        self._GM_ = GM
        self._freeze_self_()

    @lru_cache(maxsize=64)
    def elements_and_local_indices_of_dof(self, i):
        """Find the local mesh elements and local indices of global dof #i.

        When we find nothing, we return None, Otherwise, we return two lists of mesh elements, and
        local indices respectively.

        """
        mesh_elements = list()
        local_indices = list()

        for e in self._GM_: # go through all local mesh elements
            gv = self._GM_[e] # get the local gathering vector in each local mesh element
            if i in gv:
                mesh_elements.append(e)
                local_indices.append(gv.index(i))

        if len(mesh_elements) == 0:
            return None
        else:
            return mesh_elements, local_indices