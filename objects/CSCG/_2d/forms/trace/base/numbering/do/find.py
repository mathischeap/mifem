# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/19 7:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class _2dCSCG_Trace_Numbering_DO_FIND(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._freeze_self_()




    def local_dofs_on_element_edge(self, edge):
        """"""
        t = self._numbering_._tf_

        if t.__class__.__name__ == '_2dCSCG_1Trace_Outer':
            return self.___PRIVATE_trace_1_outer___(edge)

    def ___PRIVATE_trace_1_outer___(self, edge):
        """"""
        t = self._numbering_._tf_

        px, py = t.space.p

        if edge == 'U':
            return range(0, py)
        elif edge == 'D':
            return range(py, 2*py)
        elif edge == 'L':
            return range(2*py, 2*py+px)
        else:
            return range(2*py+px, 2*py+2*px)



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
