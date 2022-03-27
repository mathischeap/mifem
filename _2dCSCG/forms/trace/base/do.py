from screws.freeze.base import FrozenOnly





class _2dCSCG_Trace_DO(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._tf_.___PRIVATE_do_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def resemble(self, *args, **kwargs):
        return self._tf_.___PRIVATE_do_resemble___(*args, **kwargs)

