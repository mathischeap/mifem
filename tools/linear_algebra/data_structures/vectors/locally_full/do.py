from screws.freeze.inheriting.frozen_only import FrozenOnly


class LocallyFullVectorDo(FrozenOnly):
    """"""
    def __init__(self, LFV):
        self._v_ = LFV
        self._freeze_self_()

    def distributed_to(self, *args, **kwargs):
        self._v_.___PRIVATE_be_distributed_to___(*args, **kwargs)