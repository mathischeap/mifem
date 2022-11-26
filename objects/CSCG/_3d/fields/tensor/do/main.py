


from components.freeze.main import FrozenOnly
from objects.CSCG._3d.fields.tensor.do.reconstruct.main import _3dCSCG_Tensor_Do_Reconstruct




class _3dCSCG_TensorField_DO(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._reconstruct_ = _3dCSCG_Tensor_Do_Reconstruct(tf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._tf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_