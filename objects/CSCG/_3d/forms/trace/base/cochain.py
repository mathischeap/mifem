
from objects.CSCG.base.forms.trace.cochain import CSCG_Trace_Form_Cochain_BASE


class _3dCSCG_Trace_Cochain(CSCG_Trace_Form_Cochain_BASE):
    def __init__(self, tf):
        super().__init__(tf)

    def ___local_2_local_TEW___(self):
        """"""
        BO = self._tf_.num.basis_onside
        INDICES = [0,]
        for sn in 'NSWEBF':
            # noinspection PyUnresolvedReferences
            INDICES.append(INDICES[-1]+BO[sn])
        _D_ = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}
        TEW = dict()
        for i in self._tf_.mesh.trace.elements:
            TE = self._tf_.mesh.trace.elements[i]
            CE, CS = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_side
            i0 = _D_[CS]
            TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
        self._local_TEW_ = TEW