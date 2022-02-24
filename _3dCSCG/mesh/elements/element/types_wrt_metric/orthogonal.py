


from SCREWS.decorators import accepts
from _3dCSCG.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase


class OrthogonalElement(ElementTypeWr2MetricBase):
    """
    For orthogonal element, we only need to know its lengths on three directions, then we can fix its metric.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(3)'))
    def __init__(self, LxLyLz):
        Lx, Ly, Lz = LxLyLz
        self._mark_ = 'Orth.' + 'x{}y{}z{}'.format('%.6f' % Lx, '%.6f' % Ly, '%.6f' % Lz)
        self._freeze_self_()
