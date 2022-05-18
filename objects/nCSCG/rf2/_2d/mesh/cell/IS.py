# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class _2nCSCG_CellIS(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._atb_ = None
        self._U = None
        self._D = None
        self._L = None
        self._R = None
        self._freeze_self_()

    @property
    def root(self):
        return self._cell_.___isroot___

    @property
    def attached_to_Lv0cell_boundary(self):
        """If this cell is attached to the cscg-element(level-0-cell)-boundary (lv0-trace-element)."""
        if self._atb_ is None:
            indices = self._cell_.indices
            LEN = len(indices)
            if LEN <= 2:
                self._atb_ = True
                if LEN == 1:
                    self._U = self._D = self._L = self._R = True
                else:
                    i1 = indices[1]
                    if i1 == 0:
                        self._U = self._L = True
                        self._D = self._R = False
                    elif i1 == 1:
                        self._D = self._L = True
                        self._U = self._R = False
                    elif i1 == 2:
                        self._U = self._R = True
                        self._D = self._L = False
                    elif i1 == 3:
                        self._D = self._R = True
                        self._U = self._L = False
                    else:
                        raise Exception()
            else:
                origin, delta = self._cell_.mesh.do.find.origin_and_delta(indices)
                end = [_ + delta for _ in origin]
                if -1 in origin:
                    self._atb_ = True
                else:
                    if 1 in end:
                        self._atb_ = True
                    else:
                        self._atb_ = False

                if origin[0] == -1:
                    self._U = True
                else:
                    self._U = False

                if end[0] == 1:
                    self._D = True
                else:
                    self._D = False

                if origin[1] == -1:
                    self._L = True
                else:
                    self._L = False

                if end[1] == 1:
                    self._R = True
                else:
                    self._R = False

        return self._atb_

    @property
    def attached_to_Lv0cell_U_boundary(self):
        if self._U is None:
            _ = self.attached_to_Lv0cell_boundary
        return self._U

    @property
    def attached_to_Lv0cell_D_boundary(self):
        if self._D is None:
            _ = self.attached_to_Lv0cell_boundary
        return self._D

    @property
    def attached_to_Lv0cell_L_boundary(self):
        if self._L is None:
            _ = self.attached_to_Lv0cell_boundary
        return self._L

    @property
    def attached_to_Lv0cell_R_boundary(self):
        if self._R is None:
            _ = self.attached_to_Lv0cell_boundary
        return self._R







if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/cell/IS.py
    pass
