
from screws.freeze.main import FrozenOnly

class ___LinearSystem_Condition___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()
