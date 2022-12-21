# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/21 11:42
"""
from components.freeze.main import FrozenOnly

import numpy as np


class ___PRIVATE_eliminate_caller___(FrozenOnly):
    """"""
    def __init__(self, MDM, dim, form):
        """"""
        self._mdm_ = MDM
        self._dim_ = dim
        self._form_ = form

        _indices = 'abcdefghijklmnopq'[:MDM.ndim]
        _eli_ind = _indices[dim]
        _eli_indices = _indices.replace(_eli_ind, '')

        self._ein_str_ = _indices+','+_eli_ind+'->'+_eli_indices

        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        local_cochain = self._form_.cochain.local[basic_unit]
        mdm_minus1 = np.einsum(self._ein_str_,
                               self._mdm_[basic_unit], local_cochain, optimize='optimal')
        return mdm_minus1
