"""
BiCGSTAB

Not working.

"""



from tools.linear_algebra.preconditioners.allocator import PreconditionerAllocator
from screws.miscellaneous.timer import MyTimer
from tools.linear_algebra.solvers.parallel.BiCGSTAB.components.mpi_v0 import ___mpi_v0_BiCGSTAB___

from tools.linear_algebra.solvers.parallel.base import ParallelSolverBase

class BiCGSTAB(ParallelSolverBase):
    """"""
    def __init__(self, routine='auto', name=None):
        """"""
        super().__init__(routine, name)

    def __call__(self, A, b, x0,
                 maxiter=20, tol=1e-3, atol=1e-4,
                 preconditioner=(None, dict()),
                 COD=True
        ):
        """
        :param A: GlobalMatrix
        :param b: GlobalVector
        :param x0: LocallyFullVector
        :param maxiter:
        :param tol: tolerance.
        :param atol: absolute tolerance.
        :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
        :param COD: Clear Original Data?
        :return: Return a tuple of 5 outputs:

                1. (DistributedVector) results -- The result vector.
                2. (int) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence


                3. (float) beta -- The residual.
                4. (int) ITER -- The number of outer iterations.
                5. (str) message

        """

        message = "BiCGSTAB-" + MyTimer.current_time()

        assert x0.__class__.__name__ == "LocallyFullVector", \
                         f"x0 needs to be a 'LocallyFullVector'. Now I get {b.__class__}."

        assert maxiter >= 1 and maxiter % 1 == 0, f"maxiter={maxiter} must be >= 1."
        assert tol > 0 and atol > 0, f"tol={tol} and atol={atol} wrong, they must be > 0."

        # -------  Decide preconditioner -----------------------------------------------------------------------------------
        if preconditioner is None: preconditioner = (None, dict())

        preconditioner_ID, preconditioner_kwargs = preconditioner
        if preconditioner_ID is not None:
            assert preconditioner_ID in PreconditionerAllocator.___defined_preconditioners___(), \
                f"preconditioner={preconditioner_ID} is not coded, try one of " \
                f"{PreconditionerAllocator.___defined_preconditioners___().keys()}"
            preconditioner = PreconditionerAllocator(preconditioner_ID)(A, **preconditioner_kwargs)
        else:
            preconditioner = None

        # -------  Decide routine ------------------------------------------------------------------------------------------
        if self._routine_ == 'auto':
            ROUTINE = ___mpi_v0_BiCGSTAB___ # in the future, we may want to make de function to decide which one is the best for particular matrices.
        else:
            if self._routine_ == '0':
                ROUTINE = ___mpi_v0_BiCGSTAB___
            else:
                raise Exception(f"routine={self._routine_} is wrong.")

        # ---------- Do the computation ------------------------------------------------------------------------------------
        results, info, beta, ITER, solver_message = \
        ROUTINE(A, b, x0,
                       maxiter=maxiter, tol=tol, atol=atol,
                       preconditioner=preconditioner,
                       COD=COD
                       )

        MESSAGE =  message + '-' + solver_message
        #===================================================================================================================

        return results, info, beta, ITER, MESSAGE



