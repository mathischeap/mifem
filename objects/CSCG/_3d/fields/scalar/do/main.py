# -*- coding: utf-8 -*-

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.fields.scalar.do.reconstruct.main import _3dCSCG_Scalar_Do_Reconstruct


class _3dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._reconstruct_ = _3dCSCG_Scalar_Do_Reconstruct(sf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_