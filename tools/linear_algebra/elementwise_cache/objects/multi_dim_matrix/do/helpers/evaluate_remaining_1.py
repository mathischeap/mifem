# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/8 20:14
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

from numpy import einsum, newaxis
from scipy.sparse import csc_matrix


class nLS_DoEvaRM1(FrozenOnly):
    def __init__(self, MDM, eliminate_dimensions):
        """"""
        self._mdm_ = MDM
        _indices = 'abcdefghijklmnopq'[:MDM.ndim]
        for i in range(MDM.ndim):
            if i not in eliminate_dimensions:
                break
            else:
                pass
        # noinspection PyUnboundLocalVariable
        _rm_dim = _indices[i]
        _eli_indices = list()
        self._eli_form = list()
        for i in eliminate_dimensions:
            _eli_indices.append(_indices[i])
            self._eli_form.append(eliminate_dimensions[i])
        self._ein_str = _indices + ',' + ','.join(_eli_indices) + '->' + _rm_dim
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        vec = einsum(self._ein_str,
               self._mdm_[basic_unit],
               *[_.cochain.local[basic_unit] for _ in self._eli_form],
               optimize='optimal')

        return csc_matrix(vec[:, newaxis])


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
