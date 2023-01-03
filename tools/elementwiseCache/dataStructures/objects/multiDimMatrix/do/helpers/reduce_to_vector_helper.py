# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/11/2022 11:38 PM
"""
from components.freeze.main import FrozenOnly
from numpy import einsum, newaxis
from scipy.sparse import csc_matrix


class nLS_DoR2V_Helper(FrozenOnly):
    def __init__(self, MDM, v):
        """"""
        self._mdm_ = MDM
        self._v_ = v
        indices = 'abcdefghijklmnopq'[:MDM.ndim]

        rm_ind = MDM.correspondence.index(v)
        rm_str = indices[rm_ind]

        eli_ids = list()
        eli_fms = list()
        for i, form in enumerate(MDM.correspondence):
            if i != rm_ind:
                eli_ids.append(indices[i])
                eli_fms.append(form)
        self._ein_str = indices + ',' + ','.join(eli_ids) + '->' + rm_str
        self._eli_fms = eli_fms
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        vec = einsum(
            self._ein_str,
            self._mdm_[basic_unit],
            *[_.cochain.local[basic_unit] for _ in self._eli_fms],
            optimize='optimal'
        )

        return csc_matrix(vec[:, newaxis])
