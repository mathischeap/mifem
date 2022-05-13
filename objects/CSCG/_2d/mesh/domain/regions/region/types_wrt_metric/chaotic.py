# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2Metric
from root.config.main import rAnk, cOmm
import numpy as np
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement

class Chaotic(TypeWr2Metric):
    """
    Chaotic regions is the default regions type. If we do not mention the regions type in the domain input file,
    we will use this regions type as its type.

    If a regions is classified as a chaotic regions, then all the elements in this regions will also be chaotic.
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
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(2)]), f"spacing={spacing} is wrong."
        return ChaoticElement()


