




from screws.decorators.accepts import accepts
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase


class OrthogonalElement(ElementTypeWr2MetricBase):
    """
    An orthogonal element is an element: 1) it is a rectangle (including square), 2) its internal
    transformation is linear, 3) its left edge is parallel with x-axis.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(2)'))
    def __init__(self, LxLy):
        Lx, Ly = LxLy
        self._mark_ = 'Orth.' + 'x{}y{}'.format('%.4f' % Lx, '%.4f' % Ly)
        self._freeze_self_()