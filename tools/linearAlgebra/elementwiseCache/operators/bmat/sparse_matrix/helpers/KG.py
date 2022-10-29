
class ___BMAT_HELPER_KeyGenerator___:
    """"""
    def __init__(self, blocks):
        self.valid_block = list()
        self._check_repeated_ = True
        self._IS_CT_ = False
        self._IS_NC_ = False
        self._CT_ = '>CT<'
        self._NC_ = '>NC<'
        for i, Bi in enumerate(blocks):
            for j, Bij in enumerate(Bi):
                if Bij is not None:
                    self.valid_block.append(Bij)

    def __call__(self, item):

        if self._check_repeated_:
            self._check_repeated_ = False

            key_out = ''
            for vb in self.valid_block:
                key_out += vb._KG_(item)

            if self._NC_ in key_out:
                self._IS_NC_ = True

            else:
                temp = (key_out + key_out).find(key_out, 1, -1)
                if temp != -1:
                    repeated = key_out[:temp]
                else:
                    repeated = ''

                if repeated == self._CT_:
                    self._IS_CT_ = True
                else:
                    pass
        else:
            pass

        if self._IS_CT_:
            return self._CT_

        elif self._IS_NC_:
            return self._NC_

        else:
            key_out = ''
            for vb in self.valid_block:
                key_out += vb._KG_(item)

        return key_out
