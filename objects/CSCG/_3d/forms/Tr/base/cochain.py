
from objects.CSCG.base.forms.Tr.cochain import CSCG_TrForm_Cochain_BASE


class _3dCSCG_Tr_Cochain(CSCG_TrForm_Cochain_BASE):
    def __init__(self, Tr):
        super(_3dCSCG_Tr_Cochain, self).__init__(Tr)


    def ___local_2_local_TEW___(self):
        """"""
        if self._Tr_.k == 2:
            BO = self._Tr_.num.basis_onside
            INDICES = [0,]
            for sn in 'NSWEBF':
                # noinspection PyUnresolvedReferences
                INDICES.append(INDICES[-1]+BO[sn])
            _D_ = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}
            TEW = dict()
            for i in self._Tr_.mesh.trace.elements:
                TE = self._Tr_.mesh.trace.elements[i]
                CE, CS = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_side
                i0 = _D_[CS]
                TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
            self._local_TEW_ = TEW

        else:
            raise NotImplementedError()