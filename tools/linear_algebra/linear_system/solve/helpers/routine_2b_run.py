from screws.freeze.base import FrozenOnly



class RoutineToBeRun(FrozenOnly):
    """This actually is a wrapper of the solver A and b."""
    def __init__(self, routine, A, b):
        """"""
        self._r_ = routine
        self._A_ = A
        self._b_ = b
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """Here we actually run a routine, for example see GMRES."""
        return self._r_(self._A_, self._b_, *args, **kwargs)