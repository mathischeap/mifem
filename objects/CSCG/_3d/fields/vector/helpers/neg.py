

class ___VECTOR_NEG_HELPER_1___(object):
    def __init__(self, v):
        self._v_ = v

    def __call__(self, t, x, y, z):
        return - self._v_(t, x, y, z)
