# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2Metric
import numpy as np

from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.orthogonal import OrthogonalElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.orthogonal import OrthogonalTraceElement


class Orthogonal(TypeWr2Metric):
    """"""
    def __init__(self, region):
        """"""
        super().__init__(region)

        cc = region.corner_coordinates

        NWB, SWB, NEB, SEB, NWF, SWF, NEF, SEF = cc

        NWB_x, NWB_y, NWB_z = NWB
        SWB_x, SWB_y, SWB_z = SWB
        NEB_x, NEB_y, NEB_z = NEB
        SEB_x, SEB_y, SEB_z = SEB
        NWF_x, NWF_y, NWF_z = NWF
        SWF_x, SWF_y, SWF_z = SWF
        NEF_x, NEF_y, NEF_z = NEF
        SEF_x, SEF_y, SEF_z = SEF

        if NWB_x == NEB_x == NWF_x == NEF_x and \
           NWB_y == SWB_y == NWF_y == SWF_y and \
           NWB_z == SWB_z == NEB_z == SEB_z and \
           SWB_x == SEB_x == SWF_x == SEF_x and \
           NEB_y == SEB_y == NEF_y == SEF_y and \
           NWF_z == SWF_z == NEF_z == SEF_z:

            _Length_ = SWB_x - NWB_x
            _width_ = NEB_y - NWB_y
            _height_ = NWF_z - NWB_z
            assert _Length_ > 0 and _width_ > 0 and _height_ > 0, \
                f"orthogonal regions {region.name} " \
                f"L, W, H = {_Length_}, {_width_}, {_height_} wrong!"
            self._Lxyz_ = _Length_, _width_, _height_

        else:
            raise Exception()

        self._freeze_self_()

    @property
    def mark(self):
        if self._mark_ is None:
            self._mark_ = 'orthogonal:L{}_W{}_H{}'.format(
                '%.8f' % self._Lxyz_[0], '%.8f' % self._Lxyz_[1], '%.8f' % self._Lxyz_[2])

        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple):
        assert np.shape(spacing) == (3, 2), \
            "I need a spacing of shape (3,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(3)]), \
            f"spacing={spacing} is wrong."

        LxLyLz = [(spacing[i][1] - spacing[i][0]) * self._Lxyz_[i] for i in range(3)]
        return OrthogonalElement(LxLyLz)

    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple):
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

        return OrthogonalTraceElement(perp_to, d1, d2)
