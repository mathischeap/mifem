# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 1:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.mesh.segments.segment.types_wrt_metric.base import mpRfT2_SegmentTypeWr2Metric_Base


class mpRfT2_StraightSegment(mpRfT2_SegmentTypeWr2Metric_Base):
    """
    A straight segment.
    """
    def __init__(self, angle, length):

        if isinstance(angle, int):
            angle = str(angle)
        else:
            angle = '%.4f' % angle

        if isinstance(length, int):
            length = str(length)
        else:
            length = '%.4f' % length

        self._mark_ = 'Straight.a{}l{}'.format(angle, length)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
