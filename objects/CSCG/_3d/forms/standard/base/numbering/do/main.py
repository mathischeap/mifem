# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard.base.numbering.do.find import \
    _3dCSCG_Standard_Form_Numbering_DO_FIND_


class _3dCSCG_Standard_Form_Numbering_DO_(FrozenOnly):
    """"""
    def __init__(self, numbering):
        self._numbering_ = numbering
        self._find_ = _3dCSCG_Standard_Form_Numbering_DO_FIND_(self)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_

    def reset_cache(self):
        self._numbering_.___PRIVATE_reset_cache___()