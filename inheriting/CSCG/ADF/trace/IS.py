from screws.freeze.inheriting.frozen_only import FrozenOnly


class CSCG_ADT_TF_IS(FrozenOnly):
    """"""
    def __init__(self, adt):
        self._adt_ = adt
        self._freeze_self_()

    @property
    def hybrid(self):
        return True

    @property
    def inner_oriented(self):
        return True if self._adt_._orientation_ == 'inner' else False