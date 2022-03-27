




from objects.CSCG._2d.mesh.trace.elements.element.types_wrt_metric.base import TraceElementTypeWr2MetricBase



class ChaoticTraceElement(TraceElementTypeWr2MetricBase):
    """
    Chaotic trace element is the trace element that its metric is unique. So there is not any other
    trace elements can be the same with it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()