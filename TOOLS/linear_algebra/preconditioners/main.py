

"""
A main body of all preconditioners.

"""

from SCREWS.frozen import FrozenOnly



class Preconditioner(FrozenOnly):
    """"""
    def __init__(self, A):
        assert str(A.__class__) == "<class 'TOOLS.linear_algebra.data_structures.GlobalMatrix'>", \
            f"I need a 'TOOLS.linear_algebra.data_structures.GlobalMatrix'. Now I get {A.__class__}."
        self._A_ = A

    @property
    def M(self):
        """The matrix that approximate the inverse of A (Ax = b). So we need to do M dot A, M dot b."""
        raise NotImplementedError()

    @property
    def applying_method(self):
        """An indicator that explains how to apply this preconditioner."""
        AM = self.___applying_method___

        assert AM in ('left_multiply',)

        return AM

    @property
    def ___applying_method___(self):
        """An indicator that explains how to apply this preconditioner."""
        raise NotImplementedError()