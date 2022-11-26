

from components.freeze.base import FrozenOnly



class _3dCSCG_0TF_DOF_Matplot(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""