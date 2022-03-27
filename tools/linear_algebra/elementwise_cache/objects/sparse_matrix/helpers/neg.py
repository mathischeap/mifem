



from screws.freeze.base import FrozenOnly




class ___NEG___(FrozenOnly):
    def __init__(self, ewc):
        self._ewc_ = ewc
        self._freeze_self_()

    def __call__(self, item):
        return - self._ewc_[item]