
from screws.freeze.base import FrozenOnly



class GPD_2SF(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""