# -*- coding: utf-8 -*-
from components.miscellaneous.timer import MyTimer

from tools.linearAlgebra.solvers.Schur.base import SchurSolverBase

from tools.linearAlgebra.solvers.Schur.rank2.helpers.scipy_sparse_linalg_direct import \
    ___scipy_sparse_linalg_direct___


class Rank2(SchurSolverBase):
    """
    When rank = 2, we consider the system as

    | M B | | g |
    | C D | | h |

    And the blocks must be an integer. Let's say, blocks=i. Then the M-block is the left-upper [:i, :i]
    blocks of the `bmat`-ed `EWC_SparseMatrix` matrix.

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

                1. (dict) results -- The result vector.
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

        assert A_SHAPE[0] == A_SHAPE[1], f"block shape must be square."
        assert blocks < A_SHAPE[0], f"blocks={blocks} wrong, must be < {A_SHAPE[0]}."

        AS = A_SHAPE[0]
        if AS - blocks  == 1:
            M = A.blocks[:blocks, :blocks]
            B = A.blocks[:blocks, blocks:]
            C = A.blocks[blocks:, :blocks]
            D = A.blocks[blocks, blocks]

            g = b.blocks[:blocks]
            h = b.blocks[blocks]
        else:
            M = A.blocks[:blocks, :blocks]
            B = A.blocks[:blocks, blocks:]
            C = A.blocks[blocks:, :blocks]
            D = A.blocks[blocks:, blocks:]

            g = b.blocks[:blocks]
            h = b.blocks[blocks:]

        if routine == 'auto':
            ROUTINE = ___scipy_sparse_linalg_direct___
        else:
            raise NotImplementedError()

        # ---------- Do the computation ----------------------------------------------------------------
        if ROUTINE == ___scipy_sparse_linalg_direct___:
            GM0, GM1 = A.gathering_matrices
            GM_row = GM0.GMs[blocks:]
            GM_col = GM1.GMs[blocks:]
            results, info, beta, ITER, solver_message = ROUTINE(M, B, C, D, g, h, GM_row, GM_col)
        else:
            raise Exception()

        _ = kwargs # trivial; just leave freedom for future updates for kwargs.

        MESSAGE =  message + ':' + solver_message
        #===============================================================================================

        return results, info, beta, ITER, MESSAGE
