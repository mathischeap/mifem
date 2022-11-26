# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/18 22:22
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class MDM_IS(FrozenOnly):
    def __init__(self, MDM):
        """"""
        self._MDM_ = MDM
        self._IH_ = None
        self._freeze_self_()

    @property
    def inhomogeneous(self):
        """Return Ture if all correspondence are different."""
        if self._IH_ is None:
            cor = self._MDM_.correspondence
            IH = True
            for i, c in enumerate(cor):
                for j, r in enumerate(cor):
                    if i != j:
                        if c is r:
                            IH = False
                            break
                        else:
                            pass
                if not IH:
                    break
            self._IH_ = IH

        return self._IH_




if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/elementwise_cache/objects/multi_dim_matrix/IS.py
    pass