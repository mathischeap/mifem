# -*- coding: utf-8 -*-
from objects.CSCG._3d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2MetricBase

import numpy as np
from root.config.main import cOmm, rAnk
from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.chaotic import  ChaoticTraceElement



class Chaotic(TypeWr2MetricBase):
    """
    Chaotic regions is the default regions type. If we do not mention the regions type in the domain input file,
    we will use this regions type as its type.

    If a regions is classified as a chaotic regions, then all the elements in this region will also be chaotic.
    Therefore, we say that all elements are different. As a result, when we compute, for example, the Jacobian
    of elements, we have to do it for all elements. So, we should better avoid this.
    """
    def __init__(self, region):
        super().__init__(region)
        _ = self.mark
        self._freeze_self_()

    @property
    def mark(self):
        if self._mark_ is None:
            if rAnk == 0:
                self._mark_ = 'chaotic:' + str(id(self))
            else:
                self._mark_ = None
            self._mark_ = cOmm.bcast(self._mark_, root=0)
        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple) -> ChaoticElement:
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(3)]), f"spacing={spacing} is wrong."
        return ChaoticElement()

    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple) -> ChaoticTraceElement:
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        return ChaoticTraceElement()

