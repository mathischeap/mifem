# -*- coding: utf-8 -*-
import numpy as np

from tools.miLinearAlgebra.solvers.regular.direct.main import Direct


def ___scipy_sparse_linalg_direct___(M, B, C, D, g, h, GM_row, GM_col):
    """
    | M B | | g |
    | C D | | h |

    Parameters
    ----------
    M
    B
    C
    D
    g
    h
    GM_row
    GM_col

    Returns
    -------
    X : dict -- A dict contain the result local vector.
    info : The info which provides convergence information:
        * 0 : successful exit
        * >0 : convergence to tolerance not achieved, number of iterations
        * -1 : divergence
    beta : The residual (0).
    ITER : The number of outer iterations (0).
    message : str

    """

    inv_M = M.inv

    rA = D - C @ inv_M @ B
    rb = h - C @ inv_M @ g

    rA.gathering_matrices = (GM_row, GM_col)
    rb.gathering_matrix = GM_col

    GM = rb.gathering_matrix # cannot use GM_col because we have made a chain-GM in `rb`.

    rA = rA.assembled
    rb = rb.assembled

    # print(rA.condition.condition_number)

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