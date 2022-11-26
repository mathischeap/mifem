




from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.trace.base.numbering.do.find import _3dCSCG_Trace_Numbering_DO_FIND



class _3dCSCG_Trace_Numbering_DO(FrozenOnly):
    def __init__(self, TN):
        self._numbering_ = TN
        self._find_ = _3dCSCG_Trace_Numbering_DO_FIND(self)
        self._freeze_self_()


    @property
    def find(self):
        return self._find_




