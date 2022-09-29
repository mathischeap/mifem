# -*- coding: utf-8 -*-


from screws.decorators.accepts import accepts
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase

# from objects.mpRfT._2d.mesh.cell.types_wrt_metric.orthogonal import mpRfT2_OrthogonalCell
# from objects.mpRfT._2d.mesh.segments.segment.types_wrt_metric.straight import mpRfT2_StraightSegment

class OrthogonalElement(ElementTypeWr2MetricBase):
    """
    An orthogonal element is an element: 1) it is a rectangle (including square), 2) its internal
    transformation is linear, 3) its left edge is parallel with x-axis.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(2)'))
    def __init__(self, LxLy):
        self._LxLy_ = LxLy
        Lx, Ly = LxLy
        if Lx == Ly:
            self._mark_ = 'Orth.{:.4f}'.format(Lx)
        else:
            self._mark_ = 'Orth.x{:.4f}y{:.4f}'.format(Lx, Ly)
        self._freeze_self_()











    # def ___CLASSIFY_mpRfT2_CELL_of_origin_and_delta___(self, origin_and_delta):
    #     """"""
    #     delta = origin_and_delta[1]
    #     Lx, Ly = self._LxLy_
    #     Lx *= delta / 2
    #     Ly *= delta / 2
    #     return  mpRfT2_OrthogonalCell(Lx, Ly)
    #
    # def ___CLASSIFY_mpRfT2_segment___(self, seg):
    #     """"""
    #     direction = seg.direction
    #     if direction == 'UD':
    #         angle = 0
    #         L = self._LxLy_[0]
    #     else:
    #         angle = 90
    #         L = self._LxLy_[1]
    #
    #     rp = seg.__repr__()
    #
    #     if rp[3] == 'c':
    #         ind = rp.split(':')[-1]
    #     elif rp[3] == 't':
    #         ind =  rp.split('-')[-1]
    #     else:
    #         raise Exception()
    #
    #     LEN = len(ind) -  1
    #     length = 0.5 ** LEN * L
    #
    #     return mpRfT2_StraightSegment(angle, length)
