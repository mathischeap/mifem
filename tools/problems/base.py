"""
Not useful yet!

"""

from screws.freeze.main import FrozenOnly



class ProblemBase(FrozenOnly):
    """"""
    def __init__(self):
        self._freeze_self_()

    @property
    def ndim(self):
        raise NotImplementedError()