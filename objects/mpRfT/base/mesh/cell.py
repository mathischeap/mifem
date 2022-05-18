# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/17/2022 10:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class mpRfT_Mesh_Cell_Base(FrozenOnly):
    """"""
    def __init__(self, mesh, level, indices):
        """"""
        self._mesh_ = mesh
        assert level + 1 == len(indices), f"must be the case."
        self._level_ = level
        self._indices_ = indices
        self.___sub_cells___ = None
        self.___isroot___ = True # when it is initialized, it is a root cell of course.
        self._r = None

    def __repr__(self):
        if self._r is None:
            if self.level == 0:
                self._r = str(self.indices[0])
            else:
                self._r = str(self.indices[0])+'-'
                for i in self.indices[1:]:
                    self._r += str(i)
        return self._r

    @property
    def mesh(self):
        return self._mesh_

    @property
    def level(self):
        return self._level_

    @property
    def indices(self):
        return self._indices_

    @property
    def sub_cells(self):
        return self.___sub_cells___






if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
