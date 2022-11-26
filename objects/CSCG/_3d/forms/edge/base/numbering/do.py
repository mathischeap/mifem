# -*- coding: utf-8 -*-



from components.freeze.main import FrozenOnly




class _3dCSCG_Edge_Numbering_DO(FrozenOnly):
    def __init__(self, EN):
        self._numbering_ = EN
        self._freeze_self_()