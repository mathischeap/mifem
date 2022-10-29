# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 6:26 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from objects.mpRfT._2d.mesh.cell.types_wrt_metric.base import mpRfT2_CellTypeWr2Metric_Base


class mpRfT2_ParallelogramCell(mpRfT2_CellTypeWr2Metric_Base):
    """
    A parallelogramCell Cell is a Cell that:
        1) four edges are straight lines and form a parallelogram;
        2) internal transformation is linear,
        3) it is not an Orthogonal cell.

    So, we know, it can be a rectangle or square once its left edge is not parallel with x-axis. If
    its left edge is parallel with x-axis, such a rectangle or square is an orthogonal element.

    """
    def __init__(self, angleL, L, L_angle_U, U):
        self._data_ = angleL, L, L_angle_U, U
        assert 0 < L_angle_U < np.pi, f"L.U angle = {L_angle_U} wrong!"
        assert 0 <= angleL < 2 * np.pi, f" angle L={angleL} must be < 2pi."
        self._mark_ = 'Parallelogram.' + 'aL{}_L{}_{}_U{}'.format(
            '%.5f' % angleL, '%.5f' % L, '%.5f' % L_angle_U, '%.5f' % U)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
