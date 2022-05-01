import numpy as np

from tools.linear_algebra.solvers.regular.direct.main import Direct


def ___scipy_sparse_linalg_direct___(A, b, blocks):
    """
    | M B | | g |
    | C D | | h |

    :param A:
    :param b:
    :param blocks: A positive integer. We consider A[:blocks][:blocks] as a single local block.
    :returns: Return a tuple of 5 outputs:

                1. dict -- A dict contain the result local vector.
                2. (0) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence


                3. (0) beta -- The residual.
                4. (0) ITER -- The number of outer iterations.
                5. (str) message
    """

    M = A.blocks[:blocks, :blocks]
    B = A.blocks[:blocks, blocks:]
    C = A.blocks[blocks:, :blocks]
    D = A.blocks[blocks:, blocks:]

    g = b.blocks[:blocks]
    h = b.blocks[blocks:]

    inv_M = M.inv

    GM0, GM1 = A.gathering_matrices
    GM_row = GM0.GMs[blocks:]
    GM_col = GM1.GMs[blocks:]

    rA = D - C @ inv_M @ B
    rb = h - C @ inv_M @ g

    rA.gathering_matrices = (GM_row, GM_col)
    rb.gathering_matrix = GM_col

    GM = rb.gathering_matrix

    rA = rA.assembled
    rb = rb.assembled


    reduce_results, info, beta, ITER, solver_message = Direct()(rA, rb)

    del rA, rb

    X = dict()
    for i in GM:
        local_dofs = GM[i]

        Li = reduce_results.V[local_dofs]

        inv_Mi = inv_M[i]
        gi = g[i]
        Bi = B[i]

        xi = inv_Mi @ ( gi.toarray().ravel() - Bi @ Li)

        X[i] = np.concatenate([xi, Li])

    return X, info, beta, ITER, 'reduce-system=' + solver_message