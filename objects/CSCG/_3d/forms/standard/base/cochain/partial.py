# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly


class _3dCSCG_Standard_Form_Cochain_Partial(FrozenOnly):
    """For accessing a partial of the cochain."""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()