






from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase


class ChaoticElement(ElementTypeWr2MetricBase):
    """
    Chaotic element is the element that its metric is unique. So there is not any other elements can be the same with
    it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()