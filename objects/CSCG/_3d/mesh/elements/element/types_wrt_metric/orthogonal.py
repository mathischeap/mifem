# -*- coding: utf-8 -*-


from screws.decorators.accepts import accepts
from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase

from objects.nCSCG.rf2._3d.mesh.cell.types_wrt_metric.orthononal import _3nCSCG_OrthogonalCell

class OrthogonalElement(ElementTypeWr2MetricBase):
    """
    For orthogonal element, we only need to know its lengths on three directions, then we can fix its metric.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(3)'))
    def __init__(self, LxLyLz):
        self._LxLyLz_ = LxLyLz

        Lx, Ly, Lz = LxLyLz

        if Lx == Ly == Lz:
            self._mark_ = 'Orth.{}'.format('%.3f' % Lx)
        elif Lx == Ly:
            self._mark_ = 'Orth.xy{}z{}'.format('%.4f' % Lx, '%.4f' % Lz)
        elif Lx == Lz:
            self._mark_ = 'Orth.xz{}y{}'.format('%.4f' % Lx, '%.4f' % Ly)
        elif Ly == Lz:
            self._mark_ = 'Orth.yz{}x{}'.format('%.4f' % Ly, '%.4f' % Lx)
        else:
            self._mark_ = 'Orth.x{}y{}z{}'.format('%.3f' % Lx, '%.3f' % Ly, '%.3f' % Lz)

        self._freeze_self_()



    def ___CLASSIFY_3nCSCG_RF2_CELL_of_origin_and_delta___(self, origin_and_delta):
        """"""
        delta = origin_and_delta[1]
        Lx, Ly, Lz = self._LxLyLz_
        Lx *= delta / 2
        Ly *= delta / 2
        Lz *= delta / 2
        return  _3nCSCG_OrthogonalCell(Lx, Ly, Lz)
