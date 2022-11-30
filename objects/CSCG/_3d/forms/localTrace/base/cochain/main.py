# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:28 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG.base.forms.localTrace.cochain import CSCG_LocalTrace_CochainBase


class _3dCSCG_LocalTrace_Cochain(CSCG_LocalTrace_CochainBase):
    """"""

    def __init__(self, ltf):
        """"""
        self._components_cache_ = None
        super(_3dCSCG_LocalTrace_Cochain, self).__init__(ltf)
        self._freeze_self_()


    def ___components_of_cochain_on_element_side___(self, element, side):
        """Find the cochain referring to the `side` of element `element`.

        Parameters
        ----------
        element
        side

        Returns
        -------

        """
        if self._components_cache_ is None:
            NBO = self._ltf_.num.basis_onside
            N, S, W, E, B, F = NBO['N'], NBO['S'], NBO['W'], NBO['E'], NBO['B'], NBO['F']

            self._components_cache_ = {
                'N': (0, N),
                'S': (N, N + S),
                'W': (N + S, N + S + W),
                'E': (N + S + W, N + S + W + E),
                'B': (N + S + W + E, N + S + W + E + B),
                'F': (N + S + W + E + B, N + S + W + E + B + F),
            }

        else:
            pass

        i0, i1 = self._components_cache_[side]

        return self.local[element][i0:i1]

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
