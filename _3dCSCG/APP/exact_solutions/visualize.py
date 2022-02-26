

from screws.frozen import FrozenOnly


class ExactSolution_Visualize(FrozenOnly):
    def __init__(self, es):
        self._es_ = es
        self._freeze_self_()
