
from screws.freeze.base import FrozenOnly


class _3dCSCG_Tr_form_IS(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._freeze_self_()

    @property
    def hybrid(self):
        """Tr-form must not be hybrid."""
        return False

    @property
    def inner_oriented(self):
        return True if self._Tr_._orientation_ == 'inner' else False