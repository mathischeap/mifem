from screws.freeze.base import FrozenOnly


class DistributedVectorDo(FrozenOnly):
    """"""

    def __init__(self, DV):
        self._v_ = DV
        self._freeze_self_()

    def distributed_to(self, *args, **kwargs):
        self._v_.___PRIVATE_be_distributed_to___(*args, **kwargs)

    def gather_V_to_core(self, **kwargs):
        return self._v_.___PRIVATE_gather_V_to_core___(**kwargs)