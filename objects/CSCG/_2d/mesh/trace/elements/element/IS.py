# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly



class _2dCSCG_TraceElement_IS(FrozenOnly):
    """"""
    def __init__(self, element):
        """"""
        self._element_ = element
        self._shared_by_cores_ = None
        self._freeze_self_()



    @property
    def on_mesh_boundary(self):
        return self._element_._ondb_

    @property
    def on_periodic_boundary(self):
        return self._element_._onpb_

    @property
    def shared_by_cores(self):
        """"""
        if self._shared_by_cores_ is None:
            if self.on_mesh_boundary:
                self._shared_by_cores_ = False
            else:
                if int(self._element_._p1_[:-1]) in self._element_._elements_._mesh_.elements and \
                    int(self._element_._p2_[:-1]) in self._element_._elements_._mesh_.elements:
                    self._shared_by_cores_ = False
                else:
                    self._shared_by_cores_ = True
        return self._shared_by_cores_