# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 1:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._3d.mesh.cell.types_wrt_metric.base import _3nCSCG_CellTypeWr2Metric_Base

class _3nCSCG_ChaoticCell(_3nCSCG_CellTypeWr2Metric_Base):
    """
    Chaotic cell is a cell that its metric is unique. So there is not any other cells can be the same with
    it. To make sure that happens, we use its unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
