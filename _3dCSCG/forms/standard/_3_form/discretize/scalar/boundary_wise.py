from screws.freeze.inheriting.frozen_only import FrozenOnly


class _3dCSCG_Discretize_BoundaryWise(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()
