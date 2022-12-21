# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:28 PM
"""
from objects.CSCG.base.forms.localTrace.cochain.main import CSCG_LocalTrace_CochainBase
from components.distributors import VectorDistributor


class _3dCSCG_LocalTrace_Cochain(CSCG_LocalTrace_CochainBase):
    """"""

    def __init__(self, ltf):
        """"""
        super(_3dCSCG_LocalTrace_Cochain, self).__init__(ltf)
        self._melt_self_()
        if self._ltf_.whether.hybrid:
            self._components_cache_ = None
        else:
            self._distributor_ = VectorDistributor(self._ltf_.numbering.local_gathering)
        self._freeze_self_()

    def ___local_2_local_ESW___(self):
        """"""
        if self._ltf_.whether.hybrid:  # hybrid trace forms
            BO = self._ltf_.num.basis_onside
            INDICES = [0, ]
            sns = 'NSWEBF'
            for sn in sns:
                # noinspection PyUnresolvedReferences
                INDICES.append(INDICES[-1]+BO[sn])
            _D_ = {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}
            ESW: dict[int] = dict()

            for i in self.local:
                local_i = self.local[i]
                esw_i = dict()

                for side in sns:
                    i0 = _D_[side]
                    esw_i[side] = local_i[INDICES[i0]:INDICES[i0 + 1]]
                else:
                    pass

                ESW[i] = esw_i

            self._local_ESW_ = ESW

        else:  # non-hybrid trace forms, we use distributors.
            ESW: dict[int] = dict()
            for i in self.local:
                ESW[i] = self._distributor_(self.local[i])
            self._local_ESW_ = ESW

    def ___components_of_cochain_on_element_side___(self, element, side):
        """Find the cochain referring to the `side` of element `element`.

        Parameters
        ----------
        element
        side

        Returns
        -------

        """
        if self._ltf_.whether.hybrid:
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

        else:

            if self._local_ESW_ is not None:

                return self._local_ESW_[element][side]

            else:
                local_cochain = self.local[element]

                indices = self._ltf_.numbering.local_gathering[side]

                return local_cochain[indices]
