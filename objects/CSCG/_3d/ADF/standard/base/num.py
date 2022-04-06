


from screws.freeze.base import FrozenOnly







class _3dCSCG_ADF_SF_NUM(FrozenOnly):
    """"""
    def __init__(self, adf):
        self._adf_ = adf
        self._freeze_self_()

    @property
    def GLOBAL_dofs(self):
        return self._adf_.prime.num.GLOBAL_dofs