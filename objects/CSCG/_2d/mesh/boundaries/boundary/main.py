# -*- coding: utf-8 -*-

import sys
if './' not in sys.path:
    sys.path.append('/')
from components.freeze.main import FrozenOnly


class _2dCSCG_Mesh_Boundary(FrozenOnly):
    def __init__(self, bdrs, name):
        self._bdrs_ = bdrs
        self._name_ = name
        self._freeze_self_()
