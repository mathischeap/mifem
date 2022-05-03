# -*- coding: utf-8 -*-
"""
A main body of all preconditioners.

"""

from screws.freeze.main import FrozenOnly



class Preconditioner(FrozenOnly):
    """"""
    def __init__(self, A):
        assert A.__class__.__name__ == "GlobalMatrix", f"I need a 'GlobalMatrix'. Now I get {A.__class__}."
        self._A_ = A

    @property
    def applying_method(self):
        """An indicator that explains how to apply this preconditioner."""
        AM = self.___applying_method___

        assert AM in ('left_multiply_invM', # invM @ A x = invM @ b. So self must have `invM` property
                      )

        return AM

    @property
    def ___applying_method___(self):
        """An indicator that explains how to apply this preconditioner."""
        raise NotImplementedError()