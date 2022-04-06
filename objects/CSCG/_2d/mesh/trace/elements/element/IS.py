from screws.freeze.base import FrozenOnly



class _2dCSCG_TraceElement_IS(FrozenOnly):
    """"""
    def __init__(self, element):
        """"""
        self._element_ = element
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
        if self.on_mesh_boundary:
            return False
        else:
            if int(self._element_._p1_[:-1]) in self._element_._elements_._mesh_.elements and \
                int(self._element_._p2_[:-1]) in self._element_._elements_._mesh_.elements:
                return False
            else:
                return True