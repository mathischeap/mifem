from screws.freeze.base import FrozenOnly



class CSCG_FTWBC_DO(FrozenOnly):
    def __init__(self, twBC):
        self._twBC_ = twBC
        self._freeze_self_()

    def set_BC_body_as(self, body):
        #- set the BC func body, also do some checks here ----------------------
        self._twBC_.body = body