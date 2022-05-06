# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 1:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._3d.mesh.cell.types_wrt_metric.base import _3nCSCG_CellTypeWr2Metric_Base


class _3nCSCG_OrthogonalCell(_3nCSCG_CellTypeWr2Metric_Base):
    """
    An orthogonal cell is a cell: 1) it is a rectangle (including square), 2) its internal
    transformation is linear, 3) its edges are parallel with axes.
    """
    def __init__(self, Lx, Ly, Lz):
        if Lx == Ly == Lz:
            self._mark_ = 'Orth.{}'.format('%.4f' % Lx)
        else:
            self._mark_ = 'Orth.x{}y{}z{}'.format('%.4f' % Lx, '%.4f' % Ly, '%.4f' % Lz)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
