# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly

class mpRfT2_Mesh_Cell_SubCells(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._SC_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        """"""
        if i in self._SC_:
            pass
        elif isinstance(i, int) and 0 <= i < 4:
            self._SC_[i] = self._cell_.__class__(self._cell_.mesh,
                                                 self._cell_.level+1,
                                                 self._cell_.indices + (i,))
        else:
            raise Exception(f"2d mpRfT2_Mesh only have 4 sub-cells in each cell, "
                            f"index={i} for level: >{self._cell_.level+1}< wrong.")

        return self._SC_[i]




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/cell/sub_cells/main.py
    pass