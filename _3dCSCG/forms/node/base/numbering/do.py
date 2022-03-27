from screws.freeze.base import FrozenOnly



class _3dCSCG_Node_Numbering_DO(FrozenOnly):
    def __init__(self, NN):
        self._numbering_ = NN
        self._freeze_self_()