


from scipy import sparse as spspa
from scipy.sparse import linalg as spspalinalg
from time import time
from root.config.main import cOmm, rAnk, mAster_rank
from tools.linear_algebra.data_structures.global_matrix.main import LocallyFullVector


def ___scipy_sparse_linalg_v0___(A, b, COD=True):
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

    cOmm.barrier()
    t0 = time()
    # ...
    A = A.do.gather_M_to_core(core=mAster_rank, clean_local=COD)
    b = b.do.gather_V_to_core(core=mAster_rank, clean_local=COD)

    if rAnk == mAster_rank:
        shape = A.shape[0]
        assert shape == b.shape[0], f"A:{A.shape} does not match b{b.shape}."
        x = spspalinalg.spsolve(A, b)
        x = spspa.csr_matrix(x).T
    else:
        x = None
    x = cOmm.bcast(x, root=mAster_rank)
    x = LocallyFullVector(x)
    t3 = time()
    message = f'scipy_sparse_linalg_v0_direct = [SYSTEM-SHAPE: {A.shape}] >>> costs {int((t3-t0)*100)/100}s.'

    return x, 0, 0, 0, message