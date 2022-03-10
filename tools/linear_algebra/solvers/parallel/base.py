


from screws.freeze.inheriting.frozen_only import FrozenOnly

class ParallelSolverBase(FrozenOnly):
    """A base for all parallel solvers, direct or iterative."""
    def __init__(self, routine, name):
        """

        :param routine: Which particular routine we are going to use?
        :param name: The name of this solving process.
        """
        self._routine_ = routine
        self._name_ = name
        self._freeze_self_()