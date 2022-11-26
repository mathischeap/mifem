# -*- coding: utf-8 -*-


from scipy.sparse import linalg as spspalinalg
from time import time
from root.config.main import COMM, RANK, MASTER_RANK, MPI
from tools.linearAlgebra.dataStructures.globalMatrix.main import LocallyFullVector
import numpy as np

def ___scipy_sparse_linalg_v0___(A, b, COD=None):
    """

    :param A:
    :param b:
    :param COD: clean old data?
    :return: Return a tuple of 5 outputs:

                1. (LocallyFullVector) results -- The result vector.
                2. (0) info -- The info which provides convergence information:

                    * 0 : successful exit
                    * >0 : convergence to tolerance not achieved, number of iterations
                    * -1 : divergence

                3. (0) beta -- The residual.
                4. (0) ITER -- The number of outer iterations.
                5. (str) message
    """

    COMM.barrier()
    t0 = time()
    # ...
    if COD is None: COD = True

    num_dofs = b.shape[0]

    A = A.do.gather_M_to_core(core=MASTER_RANK, clean_local=COD)
    b = b.do.gather_V_to_core(core=MASTER_RANK, clean_local=COD)

    if RANK == MASTER_RANK:
        shape = A.shape[0]
        assert shape == b.shape[0], f"A:{A.shape} does not match b{b.shape}."
        x = spspalinalg.spsolve(A, b)
    else:
        pass

    COMM.barrier()
    if RANK != MASTER_RANK:
        x = np.empty((num_dofs,), dtype=float)
    else:
        pass

    COMM.Bcast([x, MPI.FLOAT], root=MASTER_RANK)
    x = LocallyFullVector(x)
    t3 = time()
    message = f'scipy_sparse_linalg_v0_direct = [SYSTEM-SHAPE: {A.shape}] >>> costs {int((t3-t0)*100)/100}s.'

    return x, 0, 0, 0, message