


from screws.freeze.main import FrozenOnly


class ___CV_NEG___(FrozenOnly):
    def __init__(self, V):
        self._V_ = V
        self._freeze_self_()

    def __call__(self, i):
        return - self._V_[i]