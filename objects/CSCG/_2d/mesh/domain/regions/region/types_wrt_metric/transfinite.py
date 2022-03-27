

from objects.CSCG._2d.mesh.domain.regions.region.types_wrt_metric.base import TypeWr2Metric

from screws.functions._2d_space.angle import angle
from screws.functions._2d_space.distance import distance
from screws.functions._2d_space.check_if_two_lines_parallel import if_two_lines_parallel
import numpy as np
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.chaotic import ChaoticElement
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.orthogonal import OrthogonalElement
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.parallelogram import ParallelogramElement
from typing import Union
from root.config.main import rAnk, cOmm





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
        assert np.shape(spacing) == (2,2), "I need a spacing of shape (2,2) to represent an element in a regions."
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