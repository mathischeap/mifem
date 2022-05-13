# -*- coding: utf-8 -*-


from root.config.main import cOmm, rAnk

import numpy as np

from objects.CSCG._3d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2MetricBase

from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement
from objects.CSCG._3d.mesh.elements.element.types_wrt_metric.orthogonal import OrthogonalElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.chaotic import  ChaoticTraceElement
from objects.CSCG._3d.mesh.trace.elements.element.types_wrt_metric.orthogonal import OrthogonalTraceElement

from typing import Union



class Transfinite(TypeWr2MetricBase):
    """A transfinite regions."""
    def __init__(self, region):
        super().__init__(region)
        _is_6_Plane_ = True
        for sn in 'NSWEBF':
            if self._region_._side_geometries_[sn].__class__.__name__ == 'Plane':
                pass
            else:
                _is_6_Plane_ = False
                break

        if not _is_6_Plane_:
            self._IS_chaotic_ = True

        else:
            self._IS_chaotic_ = False

            cc = self._region_.corner_coordinates
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

                self._IS_something_ = 'orthogonal'
                _Length_ = SWB_x - NWB_x
                _width_  = NEB_y - NWB_y
                _height_ = NWF_z - NWB_z
                assert _Length_ > 0 and _width_ > 0 and _height_ > 0, \
                    f"orthogonal regions {region.name} " \
                    f"L, W, H = {_Length_}, {_width_}, {_height_} wrong!"
                self._Lxyz_ = _Length_, _width_, _height_

            else:
                raise NotImplementedError(f"Not classify this type of transfinite regions yet.")

        _ = self.mark
        self._freeze_self_()

    @property
    def mark(self):
        if self._mark_ is None:
            if self._IS_chaotic_:
                if rAnk == 0:
                    self._mark_ = 'chaotic:' + str(id(self))
                else:
                    self._mark_ = None
                self._mark_ = cOmm.bcast(self._mark_, root=0)
            else:
                if self._IS_something_ == 'orthogonal':
                    self._mark_ = 'orthogonal:L{}_W{}_H{}'.format(
                        '%.8f' % self._Lxyz_[0], '%.8f' % self._Lxyz_[1], '%.8f' % self._Lxyz_[2])
                else:
                    raise NotImplementedError(f"Can not handle type={self._IS_something_} yet.")

        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple) -> Union[ChaoticElement, OrthogonalElement]:
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a regions."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(3)]), f"spacing={spacing} is wrong."
        if self._IS_chaotic_:
            return ChaoticElement()
        else:
            if self._IS_something_ == 'orthogonal':
                LxLyLz = [(spacing[i][1] - spacing[i][0]) * self._Lxyz_[i] for i in range(3)]
                return OrthogonalElement(LxLyLz)
            else:
                raise Exception()


    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple) -> Union[ChaoticTraceElement, OrthogonalTraceElement]:
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_TRACE_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )