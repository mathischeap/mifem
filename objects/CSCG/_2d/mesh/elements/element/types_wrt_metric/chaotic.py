# -*- coding: utf-8 -*-

from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase

# from objects.mpRfT._2d.mesh.cell.types_wrt_metric.chaotic import mpRfT2_ChaoticCell
# from objects.mpRfT._2d.mesh.segments.segment.types_wrt_metric.chaotic import mpRfT2_ChaoticSegment



class ChaoticElement(ElementTypeWr2MetricBase):
    """
    Chaotic element is the element that its metric is unique.

    So there is not any other elements can be the same with it. To make sure that happens,
    we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()


    # def ___CLASSIFY_mpRfT2_CELL_of_origin_and_delta___(self, origin_and_delta):
    #     """"""
    #     return mpRfT2_ChaoticCell()
    #
    # def ___CLASSIFY_mpRfT2_segment___(self, seg):
    #     return mpRfT2_ChaoticSegment()