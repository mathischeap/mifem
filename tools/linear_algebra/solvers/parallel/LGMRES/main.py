
"""
Loose GMRES.

See paper:
[A. H. BAKER et. al., A TECHNIQUE FOR ACCELERATING THE CONVERGENCE OF RESTARTED GMRES, SIAM J. MATRIX ANAL. APPL., 2005]

"""

from tools.linear_algebra.preconditioners.allocator import PreconditionerAllocator
from screws.miscellaneous.timer import MyTimer
from root.config.main import rAnk, mAster_rank
from tools.linear_algebra.solvers.parallel.LGMRES.helpers.mpi_v0 import ___mpi_v0_LGMRES___

from tools.linear_algebra.solvers.parallel.base import ParallelSolverBase
from tools.linear_algebra.data_structures.vectors.locally_full.main import LocallyFullVector


class LGMRES(ParallelSolverBase):
    """"""
    def __init__(self, routine='auto', name=None):
        """"""
        super().__init__(routine, name)


    def __call__(self, A, b, x0,
                 m=100, k=10, maxiter=20, tol=1e-3, atol=1e-4,
                 preconditioner=(None, dict()),
                 COD=True,
                 plot_residuals=False,
                 **kwargs
        ):
        """

        :param A: GlobalMatrix
        :param b: GlobalVector
        :param x0: LocallyFullVector
        :param m: restart = `m` + `k`
        :param k: restart = `m` + `k`
        :param maxiter:
        :param tol: tolerance.
        :param atol: absolute tolerance.
        :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
        :param COD: Clear Original Data?
        :param plot_residuals: bool, if we plot the residuals.
        :param kwargs: possible other args for particular routine.
        :return: Return a tuple of 5 outputs:

                1. (LocallyFullVector) results -- The result vector.
                2. (int) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence


                3. (float) beta -- The residual.
                4. (int) ITER -- The number of outer iterations.
                5. (str) message

        """
        message = "LGMRES-" + MyTimer.current_time()

        # ---- parse x0 ------------------------------------------------
        if x0 == 0: # we make it an empty LocallyFullVector
            x0 = LocallyFullVector(len(b))
        else:
            pass
        assert x0.__class__.__name__ == "LocallyFullVector", \
                         f"x0 needs to be a 'LocallyFullVector'. Now I get {x0.__class__}."
        #--------------------------------------------------------------------

        assert maxiter >= 1 and maxiter % 1 == 0, f"maxiter={maxiter} must be >= 1."
        assert m >= 3 and m % 1 == 0, f"restart={m} must be >= 3."
        assert k >= 0 and k % 1 == 0, f"restart={k} must be >= 0."
        if k == 0:
            if rAnk == mAster_rank:
                print(">>> WARNING: if k = 0, LGMRES is equivalent to GMRES and is slower than GMRES. "
                      "So please use GMRES.")
        else:
            pass
        assert tol > 0 and atol > 0, f"tol={tol} and atol={atol} wrong, they must be > 0."

        # -------  Decide preconditioner ---------------------------------------------------------------
        if preconditioner is None: preconditioner = (None, dict())

        preconditioner_ID, preconditioner_kwargs = preconditioner
        if preconditioner_ID is not None:
            preconditioner = PreconditionerAllocator(preconditioner_ID)(A, **preconditioner_kwargs)
        else:
            preconditioner = None

        # -------  Decide routine ----------------------------------------------------------------------
        if self._routine_ == 'auto':
            ROUTINE = ___mpi_v0_LGMRES___ # in the future, we may want to make a function to decide which one is the best
                                          # for particular matrices.
        else:
            if self._routine_ == '0':
                ROUTINE = ___mpi_v0_LGMRES___
            else:
                raise Exception(f"routine={self._routine_} is wrong.")


        # ---------- Do the computation ----------------------------------------------------------------
        results, info, beta, ITER, solver_message = \
        ROUTINE(A, b, x0,
                       m=m, k=k, maxiter=maxiter, tol=tol, atol=atol,
                       preconditioner=preconditioner,
                       COD=COD,
                       name=self._name_,
                       plot_residuals=plot_residuals
                    )

        _ = kwargs # trivial; just leave freedom for future updates for kwargs.

        MESSAGE =  message + '-' + solver_message
        #===============================================================================================

        return results, info, beta, ITER, MESSAGE
