



from screws.freeze.base import FrozenOnly




class ___VECMUL___(FrozenOnly):
    def __init__(self, EWC_S, EWC_V):
        self._ewc_S_ = EWC_S
        self._ewc_V_ = EWC_V
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc_S_[item] @ self._ewc_V_[item]

    def __KG_call__(self, item):
        return self._ewc_S_._KG_(item) + self._ewc_V_._KG_(item)