

import sys
if './' not in sys.path: sys.path.append('./')
from screws.frozen import FrozenOnly



class _2dCSCG_Mesh_Boundary(FrozenOnly):
    def __init__(self, bdrs, name):
        self._bdrs_ = bdrs
        self._name_ = name
        self._freeze_self_()
