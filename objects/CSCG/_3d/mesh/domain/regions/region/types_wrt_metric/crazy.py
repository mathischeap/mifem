# -*- coding: utf-8 -*-

import numpy as np
from objects.CSCG._3d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2MetricBase

from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement
from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.orthogonal import OrthogonalElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.chaotic import  ChaoticTraceElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.orthogonal import OrthogonalTraceElement
from typing import Union



class Crazy(TypeWr2MetricBase):
    """The crazy regions is the regions the crazy mesh uses."""
    def __init__(self, region):
        super().__init__(region)
        self._c_ = self._region_._domain_input_.c
        bounds = self._region_._domain_input_.bounds
        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        z0, z1 = bounds[2]
        self._Lxyz_ = (x1-x0, y1-y0, z1-z0)
        self._freeze_self_()

    @property
    def mark(self):
        if self._mark_ is None:
            self._mark_ = 'crazy:Lx{}_Ly{}_Lz{}~c{}'.format(
                '%.8f' % self._Lxyz_[0], '%.8f' % self._Lxyz_[1], '%.8f' % self._Lxyz_[2], '%.5f' % self._c_
            )
        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple) -> Union[ChaoticElement, OrthogonalElement]:
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(3)]), f"spacing={spacing} is wrong."
        if self._c_ == 0:
            LxLyLz = [(spacing[i][1] - spacing[i][0]) * self._Lxyz_[i] for i in range(3)]
            return OrthogonalElement(LxLyLz)
        else:
            return ChaoticElement()


    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple) -> Union[ChaoticTraceElement, OrthogonalTraceElement]:
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        assert len(trace_spacing) == 3, "It is not a trace spacing."
        s0, s1, s2 = trace_spacing
        if isinstance(s0, float):
            perp_to = 'x'
            d1, d2 = s1[1]-s1[0], s2[1]-s2[0]
        elif isinstance(s1, float):
            perp_to = 'y'
            d1, d2 = s2[1]-s2[0], s0[1]-s0[0]
        elif isinstance(s2, float):
            perp_to = 'z'
            d1, d2 = s0[1]-s0[0], s1[1]-s1[0]
        else:
            raise Exception()

        if self._c_ == 0:
            return OrthogonalTraceElement(perp_to, d1, d2)
        else:
            return ChaoticTraceElement()
