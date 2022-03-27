

import sys
if './' not in sys.path: sys.path.append('/')
from screws.freeze.main import FrozenOnly


class _2dCSCG_Mesh_Boundaries_Visualize(FrozenOnly):
    def __init__(self, bdrs):
        self._bdrs_ = bdrs
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    def matplot(self, usetex=False):
        raise NotImplementedError()