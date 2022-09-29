# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 8:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2dCSCG_Mesh_IS(FrozenOnly):
    """"""

    def __init__(self, element):
        """"""
        self._element_ = element
        self._internal_ = None
        self._orthogonal_ = None
        self._freeze_self_()

    @property
    def internal(self):
        """Return True if none of this element's four edges is on the mesh boundary."""
        if self._internal_ is None:
            positions = self._element_.position
            self._internal_ = True
            for pos in positions:
                if isinstance(pos, str):
                    self._internal_ = False
                    break
        return self._internal_

    @property
    def orthogonal(self):
        if self._orthogonal_ is None:
            eMark = self._element_.type_wrt_metric.mark
            if isinstance(eMark, str) and eMark[:4] == 'Orth':
                self._orthogonal_ = True
            else:
                self._orthogonal_ = False
        return self._orthogonal_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
