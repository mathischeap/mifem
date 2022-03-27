
from screws.freeze.base import FrozenOnly
from _2dCSCG.fields.vector.do.reconstruct.main import _2dCSCG_Vector_Do_Reconstruct


class _2dCSCG_VectorField_DO(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._reconstruct_ = _2dCSCG_Vector_Do_Reconstruct(vf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._vf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_