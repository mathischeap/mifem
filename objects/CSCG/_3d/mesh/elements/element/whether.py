# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _3dCSCG_MeshElement_Whether(FrozenOnly):
    """"""
    def __init__(self, element):
        """"""
        self._element_ = element
        self._internal_ = None
        self._orthogonal_ = None
        self._freeze_self_()

    @property
    def internal(self):
        """{bool}: If this mesh element is an internal mesh-element (no side is on the mesh-boundary).
        """
        if self._internal_ is None:
            position = self._element_.position
            is_internal = True
            for pos in position:
                if isinstance(pos, str):
                    is_internal = False
                    break
            self._internal_ = is_internal
        return self._internal_

    @property
    def orthogonal(self):
        if self._orthogonal_ is None:
            eMark = self._element_.type_wrt_metric.mark
            if isinstance(eMark, str) and eMark[:4] == 'Orth':
                self._orthogonal_ = True
            else:
                self._orthogonal_ = False

        return self._orthogonal_
