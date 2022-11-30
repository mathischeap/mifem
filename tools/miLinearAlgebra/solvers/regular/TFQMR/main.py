# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/20 8:42 PM
"""

from tools.miLinearAlgebra.preconditioners.allocator import PreconditionerAllocator
from components.miscellaneous.timer import MyTimer

from tools.miLinearAlgebra.solvers.regular.TFQMR.helpers.scipy_sparse_linalg_tfqmr import ___sp_sp_linalg_tfqmr___
from tools.miLinearAlgebra.dataStructures.vectors.locallyFull.main import LocallyFullVector

from tools.miLinearAlgebra.solvers.regular.base import ParallelSolverBase


class TFQMR(ParallelSolverBase):
    """"""
    def __init__(self, routine='auto', name=None):
        """To initialize a GMRES solver, we need to choose a particular routine and name it.

        Parameters
        ----------
        routine:
            {'auto', 'sp',}. They are:
                'auto':
                'sp': scipy sparse linalg gmres
        name
        """
        super().__init__(routine, name)


    def __call__(self, A, b, x0,
                 maxiter=20, tol=1e-5, atol=1e-4,
                 preconditioner=(None, dict()),
                 COD=True,
                 plot_residuals=False,
                 **kwargs
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
        :param plot_residuals: bool, if we plot the residuals.
        :param kwargs: possible other kwargs for particular routine.
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

        message = "GMRES-" + MyTimer.current_time()

        # ---- parse x0 ---------------------------------------------------------------------
        if x0 == 0: # we make it an empty LocallyFullVector
            x0 = LocallyFullVector(len(b))
        else:
            pass

        # ------------- check ----------------------------------------------------------------
        assert x0.__class__.__name__ == "LocallyFullVector", \
            f"x0 needs to be a 'LocallyFullVector'. Now I get {b.__class__}."

        if isinstance(maxiter, int):
            assert maxiter >= 1 and maxiter % 1 == 0, f"maxiter={maxiter} must be >= 1."
        elif isinstance(maxiter, str):
            MAXITER = int(maxiter)
            assert MAXITER >= 1 and MAXITER % 1 == 0, f"maxiter={maxiter} must be >= 1."
        else:
            raise Exception(f"maxiter={maxiter} is invalid")

        assert tol > 0 and atol > 0, f"tol={tol} and atol={atol} wrong, they must be > 0."

        # -------  Decide preconditioner -----------------------------------------------------------
        if preconditioner is None: preconditioner = (None, dict())

        preconditioner_ID, preconditioner_kwargs = preconditioner
        if preconditioner_ID is not None:
            preconditioner = PreconditionerAllocator(preconditioner_ID)(A, **preconditioner_kwargs)
        else:
            preconditioner = None

        # ------- Decide routine -------------------------------------------------------------------
        if self._routine_ == 'auto':
                ROUTINE = ___sp_sp_linalg_tfqmr___

        else:
            if self._routine_ == 'sp':
                ROUTINE = ___sp_sp_linalg_tfqmr___
            else:
                raise Exception(f"routine={self._routine_} is not implemented.")

        # noinspection PyTupleAssignmentBalance
        results, info, beta, ITER, solver_message = ROUTINE(
            A, b, x0,
            maxiter=maxiter, tol=tol, atol=atol,
                preconditioner=preconditioner,
                COD=COD,
                name=self._name_,
                plot_residuals=plot_residuals
                )
        _ = kwargs # trivial; just leave freedom for future updates for kwargs.

        MESSAGE =  message + '-' + solver_message
        #===========================================================================================

        return results, info, beta, ITER, MESSAGE

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
