
from tools.linear_algebra.solvers.parallel.base import ParallelSolverBase
from screws.miscellaneous.timer import MyTimer
from tools.linear_algebra.solvers.parallel.direct.helpers.scipy_sparse_linalg_v0 import ___scipy_sparse_linalg_v0___



class Direct(ParallelSolverBase):
    """"""
    def __init__(self, routine='auto', name=None):
        """"""
        super().__init__(routine, name)


    def __call__(self, A, b, COD=None,
                 **kwargs
        ):
        """

        :param A: GlobalMatrix
        :param b: GlobalVector
        :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
        :param COD: Clear Original Data?
        :param kwargs: possible other kwargs for particular routine.
        :return: Return a tuple of 5 outputs:

                1. (LocallyFullVector) results -- The result vector.
                2. (0) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence


                3. (0) beta -- The residual.
                4. (1) ITER -- The number of outer iterations.
                5. (str) message

        """
        message = "GMRES-" + MyTimer.current_time()

        # -------  Decide routine ----------------------------------------------------------------------
        if self._routine_ == 'auto':
            ROUTINE = ___scipy_sparse_linalg_v0___
            # in the future, we may want to make a function to decide which one is the best for particular matrices.
        else:
            if self._routine_ == '0':
                ROUTINE = ___scipy_sparse_linalg_v0___
            else:
                raise Exception(f"routine={self._routine_} is wrong.")

        # ---------- Do the computation ----------------------------------------------------------------
        results, info, beta, ITER, solver_message = ROUTINE(A, b, COD=COD)

        _ = kwargs # trivial; just leave freedom for future updates for kwargs.

        MESSAGE =  message + '-' + solver_message
        #===============================================================================================

        return results, info, beta, ITER, MESSAGE
