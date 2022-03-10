
from screws.freeze.inheriting.frozen_only import FrozenOnly


class CSCG_FTWF_DO(FrozenOnly):
    def __init__(self, twf):
        self._twf_ = twf
        self._freeze_self_()

    def set_func_body_as(self, *args, **kwargs):
        self._twf_.___DO_set_func_body_as___(*args, **kwargs)