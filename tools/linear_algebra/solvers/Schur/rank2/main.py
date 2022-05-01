
from screws.miscellaneous.timer import MyTimer

from tools.linear_algebra.solvers.Schur.base import SchurSolverBase

from tools.linear_algebra.solvers.Schur.rank2.helpers.scipy_sparse_linalg_direct import \
    ___scipy_sparse_linalg_direct___


class Rank2(SchurSolverBase):
    """
    When rank = 2, we consider the system as

    | M B | | g |
    | C D | | h |

    And the blocks must be an integer, lets say, blocks=i, The M-block is the left-upper [:i, :i]
    entries of the `bmat`ed `EWC_SparseMatrix` matrix.

    """
    def __init__(self, rank, blocks, routine, name):
        super(Rank2, self).__init__(rank, blocks, routine, name)

        assert rank == 2
        assert blocks % 1 == 0 and blocks > 0, f"blocks={blocks} is wrong, must be a positive integer."


    def __call__(self, A, b, **kwargs):
        """

        :param A: EWC_SparseMatrix
        :param b: EWC_ColumnVector
        :param kwargs: possible other kwargs for particular routine.
        :returns: Return a tuple of 5 outputs:

                1. (LocallyFullVector) results -- The result vector.
                2. (0) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence


                3. (0) beta -- The residual.
                4. (1) ITER -- The number of outer iterations.
                5. (str) message

        """
        message = "Schur-Rank2-" + MyTimer.current_time()

        blocks = self._blocks_
        routine = self._routine_

        A_SHAPE = A.blocks.shape
        assert blocks < A_SHAPE[0], f"blocks={blocks} wrong, must be < {A_SHAPE[0]}."

        if routine == 'auto':
            ROUTINE = ___scipy_sparse_linalg_direct___
        else:
            raise NotImplementedError()

        # ---------- Do the computation ----------------------------------------------------------------
        results, info, beta, ITER, solver_message = ROUTINE(A, b, blocks)

        _ = kwargs # trivial; just leave freedom for future updates for kwargs.

        MESSAGE =  message + ':' + solver_message
        #===============================================================================================

        return results, info, beta, ITER, MESSAGE
