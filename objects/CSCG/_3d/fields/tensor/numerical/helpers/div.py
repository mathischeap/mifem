
class ___TENSOR_DIVERGENCE_HELPER___(object):
    """"""
    def __init__(self, fx, fy, fz):
        self._fx_ = fx
        self._fy_ = fy
        self._fz_ = fz

    def __call__(self, t, x, y, z):
        return self._fx_(t, x, y, z) + self._fy_(t, x, y, z) + self._fz_(t, x, y, z)