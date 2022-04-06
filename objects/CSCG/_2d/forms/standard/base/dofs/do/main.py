from screws.freeze.base import FrozenOnly

from objects.CSCG._2d.forms.standard.base.dofs.do.find import _2dCSCG_SF_dofs_FIND




class _2dCSCG_SF_dofs_do(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _2dCSCG_SF_dofs_FIND(self._dofs_)
        return self._find_