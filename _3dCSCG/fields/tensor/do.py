


from screws.freeze.main import FrozenOnly




class _3dCSCG_TensorField_DO(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._tf_.___DO_evaluate_func_at_time___(time=time)

    def reconstruct(self, *args, **kwargs):
        return self._tf_.reconstruct(*args, **kwargs)