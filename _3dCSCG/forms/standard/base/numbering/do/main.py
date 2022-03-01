




from screws.frozen import FrozenOnly
from _3dCSCG.forms.standard.base.numbering.do.find import _3dCSCG_Standard_Form_Numbering_DO_FIND_


class _3dCSCG_Standard_Form_Numbering_DO_(FrozenOnly):
    def __init__(self, numbering):
        self._numbering_ = numbering
        self._find_ = _3dCSCG_Standard_Form_Numbering_DO_FIND_(self)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_
