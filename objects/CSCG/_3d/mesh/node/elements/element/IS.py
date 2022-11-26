# -*- coding: utf-8 -*-




from components.freeze.base import FrozenOnly


class _3dCSCG_NodeElement_IS(FrozenOnly):
    """"""
    def __init__(self, element):
        """"""
        self._element_ = element
        self._freeze_self_()


    @property
    def on_mesh_boundary(self):
        return True if len(self._element_.on_mesh_boundaries) > 0 else False

    @property
    def on_periodic_boundary(self):
        return self._element_._elements_._IS_on_periodic_boundary_[self._element_._i_]