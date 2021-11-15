





# noinspection PyUnresolvedReferences
class CSCG_Algebra_DUAL_Standard_Form:
    """"""

    @property
    def orientation(self):
        return self._orientation_

    @property
    def IS_hybrid(self):
        return True

    @property
    def IS_inner_oriented(self):
        return True if self._orientation_ == 'inner' else False

    @property
    def IS_volume_form(self):
        return self.ndim == self.k
