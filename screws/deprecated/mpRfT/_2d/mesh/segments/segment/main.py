# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.segments.segment.ct import mpRfT2_Segment_CT
from objects.mpRfT._2d.mesh.segments.segment.IS import mpRfT2_Segment_IS

import numpy as np





class mpRfT2_Segment(FrozenOnly):
    """"""
    def __init__(self, where, x_signature, y_signature):
        """"""
        if where.__class__.__name__ == 'mpRfT2_Mesh_Cell':
            assert where.level == 0, \
                f"I can only be built based on a basic-cell (cscg mesh element) or " \
                f"a basic-cell-trace-element."
        else:
            assert where.__class__.__name__ == 'mpRfT2_Mesh_BasicCells_TraceElement' , \
                f"I can only be built based on a basic-cell (cscg mesh element) or " \
                f"a basic-cell-trace-element."

        self._where_ = where
        self._xSignature_ = x_signature
        self._ySignature_ = y_signature
        self.___repr___ = None
        self._ct_ = None
        self._IS_ = None
        self._od = None
        self._direction_: str = '' # 'UD' or 'LR'
        self._neighbors_: list = list() # the two neighbors, could be two root-cells or
                                        # one root-cell plus a mesh boundary name.
        self._neighbors_N_: list = list()
        self._N_ = None
        self._type_wrt_metric_ = None
        self._character_ = None
        self._length_ = None
        self._freeze_self_()

    @property
    def N(self):
        """"""
        return self._N_

    @property
    def xSignature(self):
        """"""
        return self._xSignature_

    @property
    def ySignature(self):
        """"""
        return self._ySignature_

    @property
    def where(self):
        """"""
        return self._where_

    def __repr__(self):
        """"""
        if self.___repr___ is None:
            _s = 'SG-'

            if self.where.__class__.__name__ == 'mpRfT2_Mesh_Cell':
                _s += 'c' + str(self.where.indices[0])
                if isinstance(self.xSignature, str):
                    _s += ':y' + str(self.ySignature) + ':x' + self.xSignature
                    self._direction_ = 'UD'

                else:
                    _s += ':x' + str(self.xSignature) + ':y' + self.ySignature
                    self._direction_ = 'LR'

            elif self.where.__class__.__name__ == 'mpRfT2_Mesh_BasicCells_TraceElement':
                _s += 't' + str(self.where.i)
                if isinstance(self.xSignature, str):
                    _s += '-x' + self.xSignature
                    self._direction_ = 'UD'

                else:
                    _s += '-y' + self.ySignature
                    self._direction_ = 'LR'

            else:
                raise Exception()

            self.___repr___ = _s

        return self.___repr___

    @property
    def character(self):
        """The str representing the local position of a segment.

        For example,
            x : A UD-direction segment locally spans [-1,1]
            y1 : A LR-direction segment locally spans [0, 1]
            y01 : A LR-direction segment locally spans [-0.5, 0]

        """
        if self._character_ is None:
            rp = self.__repr__()
            if rp[3] == 'c':
                self._character_ = rp.split(':')[-1]
            elif rp[3] == 't':
                self._character_ = rp.split('-')[-1]
            else:
                raise Exception()
        return self._character_

    def __contains__(self, osg):
        """in"""
        return True if self.__repr__() in osg.__repr__() else False

    def __eq__(self, other):
        """=="""
        return True if self.__repr__() == other.__repr__() else False

    def __gt__(self, other):
        """>"""
        if self == other:
            return False
        else:
            if other in self:
                return True
            elif self in other:
                return False
            else:
                return Exception()

    def __lt__(self, other):
        """<"""
        if self == other:
            return False
        else:
            if self in other:
                return True
            elif other in self:
                return False
            else:
                raise Exception()

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = mpRfT2_Segment_IS(self)
        return self._IS_

    @property
    def coordinate_transformation(self):
        """"""
        if self._ct_ is None:
            self._ct_ = mpRfT2_Segment_CT(self)
        return self._ct_

    @property
    def origin_and_delta(self):
        """"""
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

            if self.where.__class__.__name__ == 'mpRfT2_Mesh_Cell':
                if isinstance(self.xSignature, str):
                    origin = (origin, self.ySignature)
                else:
                    origin = (self.xSignature, origin)

                # origin is a tuple representing (xi, et)-coordinates (UL corner) (in [-1,1] scale) of it.
                self._od = origin, delta

            elif self.where.__class__.__name__ == 'mpRfT2_Mesh_BasicCells_TraceElement':

                # origin is simply a float (in [-1,1] scale).
                self._od = origin, delta

            else:
                raise NotImplementedError()

        return self._od

    @property
    def direction(self):
        """'UP' or 'LR'.
        """
        if self._direction_ == '':
            _r  = self.__repr__()
        return self._direction_

    @property
    def neighbors(self):
        """The two basic cells attached to this segment. If a neighbor is not local, or it is the
        mesh boundary, we use str to represent it. Otherwise, we use the segment instance.

        The first entry is always the local segment.
        """
        return self._neighbors_

    @property
    def neighbors_N(self):
        """The N of the two neighbors. If the second neighbor is the mesh boundary, we return a
        list of length 1 representing the N of the first neighbor (must be local).

        Returns
        -------

        """
        return self._neighbors_N_

    @property
    def type_wrt_metric(self):
        if self._type_wrt_metric_ is None:
            element = self.neighbors[0].cscg_element
            self._type_wrt_metric_ = element.type_wrt_metric.___CLASSIFY_mpRfT2_segment___(self)
        return self._type_wrt_metric_

    @property
    def length(self):
        """"""
        if self._length_ is None:

            if self.IS.straight:
                x, y = self.coordinate_transformation.mapping(np.array([-1, 1]))
                x0, x1 = x
                y0, y1 = y
                self._length_ = ((x1-x0)**2 + (y1-y0)**2)**0.5
            else:
                raise NotImplementedError()

        return self._length_
    




if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/segments/segment/main.py

    from objects.mpRfT._2d.master import MeshGenerator
    mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)
    # from __init__ import rfT2
    # mesh = rfT2.rm(10, refinement_intensity=0.5)

    for sg in mesh.segments:
        # print(sg.type_wrt_metric.mark, sg.origin_and_delta)

        print(sg.length)
