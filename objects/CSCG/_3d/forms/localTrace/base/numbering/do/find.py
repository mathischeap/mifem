# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 9:27 AM
"""
from components.freeze.main import FrozenOnly
from numpy import arange


class _3dCSCG_LocalTrace_Numbering_DoFind(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._localSideCache_ = dict()
        self._freeze_self_()

    def local_dofs_on_element_side(self, side_name):
        """

        Parameters
        ----------
        side_name

        Returns
        -------

        """
        k = self._ltf_.k

        if self._ltf_.whether.hybrid:

            if k not in self._localSideCache_ is None: self._localSideCache_[k] = dict()

            if side_name not in self._localSideCache_[k]:
                NBO = self._ltf_.num.basis_onside
                NS = NBO['N']
                WE = NBO['W']
                BF = NBO['B']
                if   side_name == 'N': indices = [0, NS]
                elif side_name == 'S': indices = [NS, 2*NS]
                elif side_name == 'W': indices = [2*NS, 2*NS+WE]
                elif side_name == 'E': indices = [2*NS+WE, 2*NS+2*WE]
                elif side_name == 'B': indices = [2*NS+2*WE, 2*NS+2*WE+BF]
                elif side_name == 'F': indices = [2*NS+2*WE+BF, 2*NS+2*WE+2*BF]
                else: raise Exception()
                self._localSideCache_[k][side_name] = indices

            i0, i1 = self._localSideCache_[k][side_name]
            return arange(i0, i1)

        else:
            return self._ltf_.numbering.local_gathering[side_name]


    def dofs_on_element_side(self, element, side_name, GM=None):
        """

        Parameters
        ----------
        element
        side_name
        GM

        Returns
        -------

        """
        if GM is None:
            GM = self._ltf_.numbering.gathering
        else:
            pass

        k = self._ltf_.k

        if self._ltf_.whether.hybrid:
            return self.___Pr_find_ltf_dofs_on_element_side___(element, side_name, GM, k=k)
        else:
            local_indices = self._ltf_.numbering.local_gathering[side_name]
            return GM[element].full_vector[local_indices]


    def ___Pr_find_ltf_dofs_on_element_side___(self, element, side_name, GM, k):
        """

        Parameters
        ----------
        element
        side_name
        GM
        k

        Returns
        -------

        """
        assert self._ltf_.num.basis == GM.global_shape[1]

        if k not in self._localSideCache_ is None: self._localSideCache_[k] = dict()

        if side_name not in self._localSideCache_[k]:
            NBO = self._ltf_.num.basis_onside
            NS = NBO['N']
            WE = NBO['W']
            BF = NBO['B']
            if   side_name == 'N': indices = [0, NS]
            elif side_name == 'S': indices = [NS, 2*NS]
            elif side_name == 'W': indices = [2*NS, 2*NS+WE]
            elif side_name == 'E': indices = [2*NS+WE, 2*NS+2*WE]
            elif side_name == 'B': indices = [2*NS+2*WE, 2*NS+2*WE+BF]
            elif side_name == 'F': indices = [2*NS+2*WE+BF, 2*NS+2*WE+2*BF]
            else: raise Exception()
            self._localSideCache_[k][side_name] = indices

        i0, i1 = self._localSideCache_[k][side_name]

        return GM[element].full_vector[i0:i1]