

from screws.freeze.base import FrozenOnly

class _3dCSCG_NodeMesh_DoFind(FrozenOnly):
    """"""
    def __init__(self, NEs):
        self._NEs_ = NEs
        self._freeze_self_()