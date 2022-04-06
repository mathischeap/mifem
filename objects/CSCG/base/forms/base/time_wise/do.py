

from screws.freeze.base import FrozenOnly

class CSCG_Form_TimeWise_DO(FrozenOnly):
    def __init__(self, tw):
        self._tw_ = tw
        self._freeze_self_()

    def push_func_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_func_to_instant___(t=t)

    def push_BC_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_BC_to_instant___(t=t)

    def push_all_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_all_to_instant___(t=t)