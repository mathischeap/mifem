



from screws.freeze.base import FrozenOnly





class ___TRUE_DIV___(FrozenOnly):
    def __init__(self, ewc, number):
        self._ewc_ = ewc
        self._number_ = number
        self._freeze_self_()

    def __call__(self, item):
        return self._ewc_[item] / self._number_