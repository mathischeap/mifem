


from screws.freeze.base import FrozenOnly
from _2dCSCG.fields.scalar.do.reconstruct.main import _2dCSCG_Scalr_Do_Reconstruct


class _2dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._reconstruct_ = _2dCSCG_Scalr_Do_Reconstruct(sf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_