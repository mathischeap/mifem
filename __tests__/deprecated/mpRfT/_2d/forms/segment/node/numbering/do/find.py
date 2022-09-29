# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/15 8:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_NSgF_Numbering_DO_FIND(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._mesh_ = numbering._t_.mesh
        self._edge_names_ = ['U', 'D', 'L', 'R']
        self._cache0_ = dict()
        self._freeze_self_()

    def local_numbering_of_dofs_on(self, rc_cp, edge_name):
        """"""
        assert edge_name in self._edge_names_, f"edge_name={edge_name} is wrong."

        root_cell = self._mesh_[rc_cp]
        Frame = root_cell.frame
        start = 0
        for edge in Frame:

            end = start

            segments = Frame[edge]

            for seg in segments:
                N = self._numbering_._t_.N[seg]
                end += N + 1

            if edge == edge_name:
                break
            else:
                start = end

        return range(start, end)





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
