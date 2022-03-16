


from screws.freeze.inheriting.frozen_only import FrozenOnly




from scipy.sparse import linalg as spspalinalg

class ___LinearAlgebraINV___(FrozenOnly):
    def __init__(self, ewc):
        self._ewc_ = ewc
        self._freeze_self_()

    def __call__(self, item):
        return spspalinalg.inv(self._ewc_[item])