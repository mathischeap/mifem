
from screws.freeze.base import FrozenOnly


class _3dCSCG_ADFTr_form_IS(FrozenOnly):
    """"""
    def __init__(self, adfTr):
        self._adfTr_ = adfTr
        self._freeze_self_()

    @property
    def hybrid(self):
        """Tr-form must not be hybrid."""
        return False

    @property
    def inner_oriented(self):
        return True if self._adfTr_._orientation_ == 'inner' else False