# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/20 6:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_Mesh_InterfaceIS(FrozenOnly):
    """"""

    def __init__(self, IF):
        """"""
        self._IF_ = IF
        self._sbc_ = None
        self._freeze_self_()

    @property
    def on_mesh_boundary(self):
        return self._IF_._omb_

    @property
    def symmetric(self):
        return self._IF_._symmetric_

    @property
    def shared_by_cores(self):
        if self._sbc_ is None:
            if self._IF_.shared_with_core is not None:
                self._sbc_ = True
            else:
                self._sbc_ = False
        return self._sbc_




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
