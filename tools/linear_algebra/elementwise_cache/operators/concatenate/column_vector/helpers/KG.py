
class ___concatenate_HELPER_KeyGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors
        self._check_repeated_ = True
        self._IS_CT_ = False
        self._IS_NC_ = False
        self._CT_ = '>CT<'
        self._NC_ = '>NC<'

    def __call__(self, i):
        """"""
        if self._check_repeated_:
            self._check_repeated_ = False

            key_out = ''
            for vb in self.vectors:
                key_out += vb._KG_(i)

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

        if self._IS_CT_:
            return self._CT_

        elif self._IS_NC_:
            return self._NC_

        else:
            key_out = ''
            for vb in self.vectors:
                key_out += vb._KG_(i)

        return key_out