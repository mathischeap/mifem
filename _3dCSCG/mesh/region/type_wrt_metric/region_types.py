# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config import *
from SCREWS.frozen import FrozenOnly
from _3dCSCG.mesh.region.type_wrt_metric.element_types import ChaoticElement, OrthogonalElement
from _3dCSCG.mesh.region.type_wrt_metric.trace_element_types import  ChaoticTraceElement, OrthogonalTraceElement
from typing import Union




class TypeWr2Metric(FrozenOnly):
    def __init__(self, region):
        self._region_ = region
        self._mark_ = None

    def __eq__(self, other):
        return self.mark == other.mark and other.___IS_3dCSCG_TypeWr2Metric___
        # make sure we can comparing region type metric objects.

    @property
    def ___IS_3dCSCG_TypeWr2Metric___(self):
        return True

    @property
    def mark(self):
        raise NotImplementedError(
            f"Please implement property mark for TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple):
        """

        :param spacing: The spacing that representing the element.
        :return:
        """
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple):
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_TRACE_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )



class Chaotic(TypeWr2Metric):
    """
    Chaotic region is the default region type. If we do not mention the region type in the domain input file,
    we will use this region type as its type.

    If a region is classified as a chaotic region, then all the elements in this region will also be chaotic.
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
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a region."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(3)]), f"spacing={spacing} is wrong."
        return ChaoticElement()

    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple) -> ChaoticTraceElement:
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        return ChaoticTraceElement()




class Crazy(TypeWr2Metric):
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
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a region."
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




class Transfinite(TypeWr2Metric):
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
                    f"orthogonal region {region.name} " \
                    f"L, W, H = {_Length_}, {_width_}, {_height_} wrong!"
                self._Lxyz_ = _Length_, _width_, _height_

            else:
                raise NotImplementedError(f"Not classify this type of transfinite region yet.")

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
        assert np.shape(spacing) == (3,2), "I need a spacing of shape (3,2) to represent an element in a region."
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