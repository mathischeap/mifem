# -*- coding: utf-8 -*-


from screws.decorators.accepts import accepts
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase

from objects.nCSCG.rf2._2d.mesh.cell.types_wrt_metric.orthononal import _2nCSCG_OrthogonalCell

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

    def ___CLASSIFY_2nCSCG_RF2_CELL_of_origin_and_delta___(self, origin_and_delta):
        """"""
        delta = origin_and_delta[1]
        Lx, Ly = self._LxLy_
        Lx *= delta / 2
        Ly *= delta / 2
        return  _2nCSCG_OrthogonalCell(Lx, Ly)