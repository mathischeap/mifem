# -*- coding: utf-8 -*-


from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase

from objects.nCSCG.rf2._3d.mesh.cell.types_wrt_metric.chaotic import _3nCSCG_ChaoticCell

class ChaoticElement(ElementTypeWr2MetricBase):
    """
    Chaotic element is the element that its metric is unique. So no any other elements can be the same with
    it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()



    def ___CLASSIFY_3nCSCG_RF2_CELL_of_origin_and_delta___(self, origin_and_delta):
        """"""
        return _3nCSCG_ChaoticCell()