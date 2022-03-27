
from screws.freeze.base import FrozenOnly


class CSCG_trace_form_IS(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    @property
    def hybrid(self):
        """Trace form must be hybrid."""
        return True

    @property
    def inner_oriented(self):
        return True if self._tf_._orientation_ == 'inner' else False