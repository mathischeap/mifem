# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard.base.numbering.do.find import _2dCSCG_SF_numbering_do_find


class _2dCSCG_SF_numbering_do(FrozenOnly):
    """"""
    def __init__(self, numbering):
        self._numbering_ = numbering
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _2dCSCG_SF_numbering_do_find(self._numbering_)
        return self._find_