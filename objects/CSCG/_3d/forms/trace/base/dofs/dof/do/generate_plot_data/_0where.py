
from screws.freeze.base import FrozenOnly



class _3dCSCG_T0F_DOF_Where(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, density=10):
        """"""