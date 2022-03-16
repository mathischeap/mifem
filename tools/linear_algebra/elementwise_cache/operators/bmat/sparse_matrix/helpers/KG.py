
class ___BMAT_HELPER_KeyGenerator___:
    """"""
    def __init__(self, blocks):
        self.valid_block = list()
        for i, Bi in enumerate(blocks):
            for j, Bij in enumerate(Bi):
                if Bij is not None:
                    self.valid_block.append(Bij)

    def __call__(self, item):
        key_out = ''
        for vb in self.valid_block:
            key_out += vb._KG_(item)
        return key_out
