from screws.freeze.base import FrozenOnly



class CSCG_FTWBC_DO(FrozenOnly):
    def __init__(self, twBC):
        self._twBC_ = twBC
        self._freeze_self_()

    def set_BC_body_as(self, *args, **kwargs):
        self._twBC_.___DO_set_BC_body_as___(*args, **kwargs)


