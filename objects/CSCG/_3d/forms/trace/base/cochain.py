# -*- coding: utf-8 -*-
from objects.CSCG.base.forms.trace.cochain.main import CSCG_Trace_Form_Cochain_BASE
from components.distributors import VectorDistributor


class _3dCSCG_Trace_Cochain(CSCG_Trace_Form_Cochain_BASE):
    def __init__(self, tf):
        super().__init__(tf)
        if self._tf_.whether.hybrid:
            pass
        else:
            self._melt_self_()
            self._distributor_ = VectorDistributor(self._tf_.numbering.local_gathering)
            self._freeze_self_()

    def ___local_2_local_TEW___(self):
        """"""
        if self._tf_.whether.hybrid: # hybrid trace forms
            BO = self._tf_.num.basis_onside
            INDICES = [0,]
            for sn in 'NSWEBF':
                # noinspection PyUnresolvedReferences
                INDICES.append(INDICES[-1]+BO[sn])
            _D_ = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}
            TEW : dict[int] = dict()
            for i in self._tf_.mesh.trace.elements:
                TE = self._tf_.mesh.trace.elements[i]
                CE, CS = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_side
                if CE in self.local and i not in TEW:
                    i0 = _D_[CS]
                    TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
                else:
                    pass
            self._local_TEW_ = TEW

        else: # non-hybrid trace forms, we use distributors.
            Tmap = self._tf_.mesh.trace.elements.map
            TEW : dict[int] = dict()
            sns = 'NSWEBF'
            for ele in self.local:
                ESW = self._distributor_(self.local[ele])
                TES = Tmap[ele]
                for i, te in enumerate(TES):
                    if te in TEW:
                        continue
                    else:
                        TEW[te] = ESW[sns[i]]

            # note that it is possible that local is not full but local_TEW is full."
            self._local_TEW_ = TEW