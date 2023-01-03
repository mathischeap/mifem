

class ___TENSOR_SUB_HELPER_1___(object):
    def __init__(self, w, u):
        self._w_ = w
        self._u_ = u

    def __call__(self, t, x, y, z):
        return self._w_(t, x, y, z) - self._u_(t, x, y, z)
