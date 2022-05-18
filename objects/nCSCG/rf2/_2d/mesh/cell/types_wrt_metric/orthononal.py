# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 1:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.cell.types_wrt_metric.base import _2nCSCG_CellTypeWr2Metric_Base


class _2nCSCG_OrthogonalCell(_2nCSCG_CellTypeWr2Metric_Base):
    """
    An orthogonal cell is a cell: 1) it is a rectangle (including square), 2) its internal
    transformation is linear, 3) its edges are parallel with axes.
    """
    def __init__(self, Lx, Ly):

        Lx = '%.4f' % Lx
        Ly = '%.4f' % Ly

        if Lx == Ly:
            self._mark_ = 'Orth.{}'.format(Lx)
        else:
            self._mark_ = 'Orth.x{}y{}'.format(Lx, Ly)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
