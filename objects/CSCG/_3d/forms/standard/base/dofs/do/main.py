from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.standard.base.dofs.do.find import _3dCSCG_SF_dofs_FIND


class _3dCSCG_SF_dofs_DO(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        """"""
        if self._find_ is None:
            self._find_ = _3dCSCG_SF_dofs_FIND(self._dofs_)
        return self._find_
