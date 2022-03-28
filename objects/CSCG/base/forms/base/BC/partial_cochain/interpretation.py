


from screws.freeze.main import FrozenOnly

class _PartialCochain_Interpretation_(FrozenOnly):
    """"""
    def __init__(self, pc):
        self._pc_ = pc
        self._freeze_self_()


    @property
    def local_cochains(self):
        """The cochain property actually is already interpreted as local cochains."""
        return self._pc_.cochain