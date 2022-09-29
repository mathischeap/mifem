# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2Metric

from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.orthogonal import OrthogonalElement

from objects.CSCG._2d.mesh.trace.elements.element.types_wrt_metric.chaotic import  ChaoticTraceElement
from objects.CSCG._2d.mesh.trace.elements.element.types_wrt_metric.orthogonal import OrthogonalTraceElement

from typing import Union
import numpy as np

class Crazy(TypeWr2Metric):
    """"""
    def __init__(self, region):
        """"""
        super().__init__(region)
        self._c_ = self._region_._domain_input_.c
        bounds = self._region_._domain_input_.bounds
        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        self._Lxy_ = (x1-x0, y1-y0)
        self._freeze_self_()

    @property
    def mark(self):
        if self._mark_ is None:
            self._mark_ = 'crazy:Lx{}_Ly{}~c{}'.format(
                '%.8f' % self._Lxy_[0], '%.8f' % self._Lxy_[1], '%.8f' % self._c_
            )
        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple) -> Union[ChaoticElement, OrthogonalElement]:
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(2)]), f"spacing={spacing} is wrong."
        if self._c_ == 0:
            LxLy = [(spacing[i][1] - spacing[i][0]) * self._Lxy_[i] for i in range(2)]
            return OrthogonalElement(LxLy)
        else:
            return ChaoticElement()


    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple) -> Union[ChaoticTraceElement, OrthogonalTraceElement]:
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        if self._c_ == 0:
            assert len(trace_spacing) == 2, "It is not a trace spacing."
            s0, s1 = trace_spacing
            if isinstance(s0, float):
                along = 'y'
                d = s1[1]-s1[0]
            elif isinstance(s1, float):
                along = 'x'
                d = s0[1]-s0[0]
            else:
                raise Exception()
            return OrthogonalTraceElement(along, d)
        else:
            return ChaoticTraceElement()