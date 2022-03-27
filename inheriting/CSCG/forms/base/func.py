from screws.freeze.base import FrozenOnly




class CSCG_Form_Func(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._body_ = None
        self._ftype_ = None

    @property
    def body(self):
        """
        The function body.
        """
        return self._body_

    @property
    def ftype(self):
        """
        The function type.

        :return: The function type.
        :rtype: str
        """
        return self._ftype_

