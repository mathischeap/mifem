"""
Not useful yet!

"""

from screws.frozen import FrozenOnly



class ProblemBase(FrozenOnly):
    """"""
    def __init__(self):
        self._freeze_self_()

    @property
    def ndim(self):
        raise NotImplementedError()







class _3d_icpsNS(ProblemBase):
    """"""
    def __init__(self):
        super(_3d_icpsNS, self).__init__()

    @property
    def ndim(self):
        return 3