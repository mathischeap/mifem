from screws.freeze.inheriting.frozen_only import FrozenOnly


class CSCG_ADF_SF_IS(FrozenOnly):
    """"""
    def __init__(self, adf):
        self._adf_ = adf
        self._freeze_self_()

    @property
    def inner_oriented(self):
        return True if self._adf_._orientation_ == 'inner' else False

    @property
    def volume_form(self):
        """"""
        return self._adf_.ndim == self._adf_.k

    @property
    def hybrid(self):
        """An AD standard form must be a hybrid form."""
        return True