# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly



class _3dCSCG_TraceElement_Whether(FrozenOnly):
    """"""
    def __init__(self, element):
        self._element_ = element
        self._sbc_ = None
        self._orthogonal_ = None
        self._freeze_self_()


    @property
    def orthogonal(self):
        if self._orthogonal_ is None:
            tMark = self._element_.type_wrt_metric.mark

            if isinstance(tMark, str) and tMark[:4] == 'Orth':
                self._orthogonal_ = True
            else:
                self._orthogonal_ = False

        return self._orthogonal_

    @property
    def on_periodic_boundary(self):
        """As this property name says."""
        return self._element_._onpb_

    @property
    def shared_by_cores(self):
        """True or False, as this property name says."""
        if self._sbc_ is None:
            if self.on_mesh_boundary:
                self._sbc_ = False
            else:
                if int(self._element_._p1_[:-1]) in self._element_._elements_._mesh_.elements and \
                    int(self._element_._p2_[:-1]) in self._element_._elements_._mesh_.elements:
                    self._sbc_ = False
                else:
                    self._sbc_ = True
        return self._sbc_

    @property
    def on_mesh_boundary(self):
        """As this property name says."""
        return self._element_._ondb_