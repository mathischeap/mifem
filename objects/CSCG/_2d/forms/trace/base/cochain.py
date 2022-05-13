# -*- coding: utf-8 -*-

from objects.CSCG.base.forms.trace.cochain import CSCG_Trace_Form_Cochain_BASE


class _2dCSCG_Trace_Cochain(CSCG_Trace_Form_Cochain_BASE):
    def __init__(self, tf):
        super().__init__(tf)

    def ___local_2_local_TEW___(self):
        """"""
        BO = self._tf_.NUM_basis_onside
        INDICES = [0,]
        for sn in 'UDLR':
            # noinspection PyUnresolvedReferences
            INDICES.append(INDICES[-1]+BO[sn])
        _D_ = {'U':0, 'D':1, 'L':2, 'R':3}
        TEW = dict()
        for i in self._tf_.mesh.trace.elements:
            TE = self._tf_.mesh.trace.elements[i]
            CE, Ce = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_edge
            i0 = _D_[Ce]
            TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
        self.local_TEW = TEW
