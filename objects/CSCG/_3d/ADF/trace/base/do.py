



from screws.freeze.main import FrozenOnly


class _3dCSCG_Algebra_DUAL_Trace_Form_DO(FrozenOnly):
    def __init__(self, dt):
        self._dt_ = dt
        self._freeze_self_()


    def reset_cache(self):
        self._dt_.___PRIVATE_reset_cache___()