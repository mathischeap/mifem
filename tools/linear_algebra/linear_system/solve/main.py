# -*- coding: utf-8 -*_


from screws.freeze.main import FrozenOnly
from tools.linear_algebra.solvers.Schur.allocator import SchurSolverDistributor
from tools.linear_algebra.solvers.regular.allocator import RegularSolverDistributor
from tools.linear_algebra.linear_system.solve.helpers.routine_2b_run import RoutineToBeRun






class ___LinearSystem_Solve___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()


    def __call__(self, solver_name, **kwargs):
        """
        Assemble and solve self.

        We first assemble self into a global system and then select the solver here! Note that
        we do not really do the solving here.

        :return: Return a tuple of 5 outputs:

                0. (LocallyFullVector) results -- The result vector.
                1. (int) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence

                2. (float) beta -- The residual.
                3. (int) ITER -- The number of outer iterations.
                4. (str) message

        """
        RSN = RegularSolverDistributor.___solver_name___()
        if solver_name in RSN:
            A = self._LS_.A.assembled
            b = self._LS_.b.assembled
            solver_2b_run_with_parameters = RegularSolverDistributor(solver_name, **kwargs)
            RTB = RoutineToBeRun(solver_2b_run_with_parameters, A, b) # wrap A, b with the solver
            # call RTB with any other args and kwargs in fact means a call of solver(A, b, *args, **kwargs)
            return RTB
        elif solver_name == 'Schur':
            assert 'rank' in kwargs, "we need rank"
            assert 'blocks' in kwargs, "we need blocks"
            rank = kwargs['rank']
            blocks = kwargs['blocks']
            del kwargs['rank'], kwargs['blocks']
            solver_2b_run_with_parameters = SchurSolverDistributor(rank, blocks, **kwargs)
            RTB = RoutineToBeRun(solver_2b_run_with_parameters, self._LS_.A, self._LS_.b)
            return RTB

        else:
            raise Exception(f"solver_name={solver_name} not found!")