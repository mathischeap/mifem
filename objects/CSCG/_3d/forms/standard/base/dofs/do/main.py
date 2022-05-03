from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.standard.base.dofs.do.find import _3dCSCG_SF_dofs_FIND

from objects.CSCG._3d.forms.standard.base.dofs.do.helpers._01f_hybrid_pairing_check_ import _01f_hybrid_pairing_check_
from objects.CSCG._3d.forms.standard.base.dofs.do.helpers._2f_hybrid_pairing_check_ import _2f_hybrid_pairing_check_


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

    def hybrid_pairing_check(self, *args, **kwargs):
        sf = self._dofs_._sf_
        k = sf.k
        if k == 3:
            raise Exception(f"3-sf does not need hybridization!")
        elif k == 2:
            return _2f_hybrid_pairing_check_(sf, *args, **kwargs)
        else:
            return _01f_hybrid_pairing_check_(sf, *args, **kwargs)