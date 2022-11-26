from components.freeze.base import FrozenOnly



class _2dCSCG_SF_DOF_VIS(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._freeze_self_()
