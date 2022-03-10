
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.solvers.parallel.allocator import ParallelSolverDistributor
from tools.linear_algebra.linear_system.solve.components.routine_2b_run import RoutineToBeRun

class ___LinearSystem_Solve___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()


    def __call__(self, solver_name, routine='auto', name=''):
        """
        Assemble and solve self.

        We first assemble self into a global system and then select the solver here! Note that
        we do not really do the solving here.

        """
        A = self._LS_.A.assembled
        b = self._LS_.b.assembled
        solver_2b_run_with_parameters = ParallelSolverDistributor(solver_name, routine=routine, name=name)
        RTB = RoutineToBeRun(solver_2b_run_with_parameters, A, b)
        return RTB