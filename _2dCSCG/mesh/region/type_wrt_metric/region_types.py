# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config import *
from SCREWS.frozen import FrozenOnly
from _2dCSCG.mesh.region.type_wrt_metric.element_types import ChaoticElement
from _2dCSCG.mesh.region.type_wrt_metric.element_types import OrthogonalElement
from _2dCSCG.mesh.region.type_wrt_metric.element_types import ParallelogramElement
from typing import Union
from SCREWS.functions._2d import distance, angle, if_two_lines_parallel



class TypeWr2Metric(FrozenOnly):
    def __init__(self, region):
        self._region_ = region
        self._mark_ = None

    def __eq__(self, other):
        return self.mark == other.mark and other.___IS_2dCSCG_TypeWr2Metric___
        # make sure we can comparing region type metric objects.

    @property
    def ___IS_2dCSCG_TypeWr2Metric___(self):
        return True

    @property
    def mark(self):
        raise NotImplementedError(
            f"Please implement property mark for TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple):
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )




# -----------------------------------------------------------------------------------------------------------




class Chaotic(TypeWr2Metric):
    """
    Chaotic region is the default region type. If we do not mention the region type in the domain input file,
    we will use this region type as its type.

    If a region is classified as a chaotic region, then all the elements in this region will also be chaotic.
    Therefore, we say that all elements are different. As a result, when we compute, for example, the Jacobian
    of elements, we have to do it for all elements. So, we should better aviod this.
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
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a region."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(2)]), f"spacing={spacing} is wrong."
        return ChaoticElement()






class Crazy(TypeWr2Metric):
    def __init__(self, region):
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
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a region."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(2)]), f"spacing={spacing} is wrong."
        if self._c_ == 0:
            LxLy = [(spacing[i][1] - spacing[i][0]) * self._Lxy_[i] for i in range(2)]
            return OrthogonalElement(LxLy)
        else:
            return ChaoticElement()



class Transfinite(TypeWr2Metric):
    def __init__(self, region):
        super().__init__(region)
        _is_4_Straight_ = True
        for en in 'UDLR':
            if self._region_._edge_geometries_[en].__class__.__name__ == 'Straight':
                pass
            else:
                _is_4_Straight_ = False
                break
        self._IS_chaotic_ = False
        if not _is_4_Straight_:
            self._IS_chaotic_ = True
        else:
            cc = self._region_.corner_coordinates
            UL, DL, UR, DR = cc
            self._len_L = distance(UL, DL)
            self._len_R = distance(UR, DR)
            self._len_U = distance(UL, UR)
            self._len_D = distance(DL, DR)
            if not (self._len_L == self._len_R and self._len_U == self._len_D):
                self._IS_chaotic_ = True
            else:
                assert if_two_lines_parallel(UL, DL, UR, DR) and \
                       if_two_lines_parallel(UL, UR, DL, DR), "no interaction!"
                self._LL_ = (self._len_L, self._len_U)
                angle_L = angle(UL, DL)
                angle_U = angle(UL, UR)

                if angle_U < angle_L: angle_U += 2*np.pi

                self.___L_angle_U___ = angle_U - angle_L
                self._angle_L = angle_L

                assert self._angle_L < 2*np.pi, f" angle L={self._angle_L} must be < 2pi."
                assert 0 < self.___L_angle_U___ < np.pi, f"L.U angle = {self.___L_angle_U___} wrong!"

                try:
                    np.testing.assert_almost_equal(self.___L_angle_U___, np.pi/2)
                    assert self._angle_L % (2*np.pi) < 1e-12
                    self._IS_orthogonal_or_parallelogram_ = 'orthogonal'
                except AssertionError:
                    self._IS_orthogonal_or_parallelogram_ = 'parallelogram'
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
                if self._IS_orthogonal_or_parallelogram_ == 'orthogonal':
                    self._mark_ = 'orthogonal:UD{}_LR{}'.format(
                        '%.8f' % self._len_U, '%.8f' % self._len_L)
                elif self._IS_orthogonal_or_parallelogram_ == 'parallelogram':
                    self._mark_ = 'parallelogram:angleL{}_lenL{}_A{}A_lenU{}'.format('%.8f' % self._angle_L,
                        '%.8f' % self._len_L, '%.8f' % self.___L_angle_U___, '%.8f' % self._len_U)
                else:
                    raise Exception("Something is wrong.")
        return self._mark_

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple) -> \
        Union[ChaoticElement, OrthogonalElement,ParallelogramElement]:
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a region."
        assert all([0 <= spacing[i][0] < spacing[i][1] <= 1 for i in range(2)]), f"spacing={spacing} is wrong."
        if self._IS_chaotic_:
            return ChaoticElement()
        elif self._IS_orthogonal_or_parallelogram_ == 'orthogonal':
            LxLy = [(spacing[i][1] - spacing[i][0]) * self._LL_[i] for i in range(2)]
            return OrthogonalElement(LxLy)
        elif self._IS_orthogonal_or_parallelogram_ == 'parallelogram':
            LR_UD = [(spacing[i][1] - spacing[i][0]) * self._LL_[i] for i in range(2)]
            return ParallelogramElement(self._angle_L, LR_UD[0], self.___L_angle_U___ ,LR_UD[1])
        else:
            raise Exception("Something is wrong.")