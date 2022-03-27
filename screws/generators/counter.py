
from screws.freeze.base import FrozenOnly

class Counter(FrozenOnly):
    """"""
    def __init__(self, start=0):
        """ """
        assert start % 1 == 0
        self._start_ = int(start)
        self._i_ = int(start - 1)
        self._freeze_self_()

    def __next__(self):
        self._i_ += 1
        return self._i_






if __name__ == '__main__':
    C = Counter()
    print(next(C))
    print(next(C))
    print(next(C))
    print(next(C))