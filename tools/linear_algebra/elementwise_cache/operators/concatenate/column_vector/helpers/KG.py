
class ___concatenate_HELPER_KeyGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors

    def __call__(self, i):
        """"""
        key_out = ''
        for v in self.vectors:
            key_out += v._KG_(i)
        return key_out