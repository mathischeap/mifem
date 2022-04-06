

from screws.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.elements.do.find import _2dCSCG_Mesh_Elements_DO_FIND



class _2dCSCG_Mesh_Elements_do(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
             self._find_ = _2dCSCG_Mesh_Elements_DO_FIND(self._elements_)
        return self._find_