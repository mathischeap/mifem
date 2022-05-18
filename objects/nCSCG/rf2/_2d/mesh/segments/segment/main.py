# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 3:38 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.exceptions import _2nCSCG_RF2_SegmentComparisonError
from objects.nCSCG.rf2._2d.mesh.segments.segment.coordinate_transformation import SegmentCT

class _2nCSCG_Segment(FrozenOnly):
    """"""

    def __init__(self, where, x_signature, y_signature):
        """"""
        if where.__class__.__name__ == '_2nCSCG_RF2_MeshCell':
            assert where.level == 0, \
                f"I can only be built based on a lv0-cell (cscg mesh element) or a Lv0-trace-element."
        else:
            assert where.__class__.__name__ == '_2nCSCG_Lv0TraceElement' , \
                f"I can only be built based on a lv0-cell (cscg mesh element) or a Lv0-trace-element."
        self._where_ = where
        self._xSignature_ = x_signature
        self._ySignature_ = y_signature
        self.___repr___ = None
        self._ct_ = None
        self._od = None
        self._freeze_self_()

    @property
    def xSignature(self):
        return self._xSignature_

    @property
    def ySignature(self):
        return self._ySignature_

    @property
    def where(self):
        return self._where_

    def __repr__(self):
        """"""
        if self.___repr___ is None:
            _s = 'SG-'

            if self.where.__class__.__name__ == '_2nCSCG_RF2_MeshCell':
                _s += 'c' + str(self.where.indices[0])
                if isinstance(self.xSignature, str):
                    _s += ':y' + str(self.ySignature) + ':x' + self.xSignature
                else:
                    _s += ':x' + str(self.xSignature) + ':y' + self.ySignature

            elif self.where.__class__.__name__ == '_2nCSCG_Lv0TraceElement':
                _s += 't' + str(self.where.i)
                if isinstance(self.xSignature, str):
                    _s += '-x' + self.xSignature
                else:
                    _s += '-y' + self.ySignature

            else:
                raise Exception()

            self.___repr___ = _s
        return self.___repr___

    def __contains__(self, osg):
        """"""
        return True if self.__repr__() in osg.__repr__() else False

    def __eq__(self, other):
        """"""
        return True if self.__repr__() == other.__repr__() else False

    def __gt__(self, other):
        """"""
        if self == other:
            return False
        else:
            if other in self:
                return True
            elif self in other:
                return False
            else:
                return _2nCSCG_RF2_SegmentComparisonError()

    def __lt__(self, other):
        if self == other:
            return False
        else:
            if self in other:
                return True
            elif other in self:
                return False
            else:
                raise _2nCSCG_RF2_SegmentComparisonError()

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = SegmentCT(self)
        return self._ct_

    @property
    def origin_and_delta(self):
        if self._od is None:
            if isinstance(self.xSignature, str):
                signature = self.xSignature
            elif isinstance(self.ySignature, str):
                signature = self.ySignature
            else:
                raise Exception()

            origin, delta = -1, 2
            for s in signature:
                delta *= 0.5
                if s == '0':
                    pass
                elif s == '1':
                    origin += delta
                else:
                    raise Exception()

            if self.where.__class__.__name__ == '_2nCSCG_RF2_MeshCell':
                if isinstance(self.xSignature, str):
                    origin = (origin, self.ySignature)
                else:
                    origin = (self.xSignature, origin)

                # origin is a tuple representing (xi, et)-coordinates (in [-1,1] scale) of it.
                self._od = origin, delta


            elif self.where.__class__.__name__ == '_2nCSCG_Lv0TraceElement':

                # origin is simply a float (in [-1,1] scale).
                self._od = origin, delta

            else:
                raise NotImplementedError()

        return self._od







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/segments/segment/main.py
    pass
