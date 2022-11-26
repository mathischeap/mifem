# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/20 8:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from root.config.main import COMM, RANK, MASTER_RANK, np, MPI
from time import time

from scipy.sparse import linalg as spspalinalg
from tools.linearAlgebra.dataStructures.globalMatrix.main import LocallyFullVector



def ___sp_sp_linalg_tfqmr___(A, b, x0,
                             maxiter=20, tol=1e-3, atol=1e-4,
                             preconditioner=None,
                            COD=True, name=None, plot_residuals=False):
    """

    Parameters
    ----------
    A
    b
    x0
    maxiter
    tol
    atol
    preconditioner
    COD

    Returns
    -------

        1. (LocallyFullVector) results -- The result vector.
        2. (int) info -- The info which provides convergence information:

            * 0 : successful exit
            * >0 : (convergence to tolerance not achieved) `Info` means number of iterations
            * <0 : illegal input or breakdown

        3. (float) beta -- The residual.
        4. (int) ITER -- The number of outer iterations.
        5. (str) message

    """

    COMM.barrier()
    t0 = time()

    num_dofs = b.shape[0]

    A = A.do.gather_M_to_core(core=MASTER_RANK, clean_local=COD)
    b = b.do.gather_V_to_core(core=MASTER_RANK, clean_local=COD)

    if RANK == MASTER_RANK:
        shape = A.shape[0]
        assert shape == b.shape[0], f"A:{A.shape} does not match b{b.shape}."
        x0 = np.array(x0.V)

        if preconditioner is not None:
            if preconditioner.__class__.__name__ == 'spiLU':

                sA_iLU = spspalinalg.spilu(A,
                                           drop_tol=preconditioner.drop_tol,
                                           fill_factor=preconditioner.fill_factor,
                                           drop_rule=preconditioner.drop_rule)

                M = spspalinalg.LinearOperator(A.shape, sA_iLU.solve)
                # noinspection PyTypeChecker
                RES = spspalinalg.tfqmr(A, b,
                                      x0=x0, tol=tol, maxiter=maxiter,
                                      M=M,
                                      atol=atol
                                      )
            else:
                raise Exception(f"Cannot used {preconditioner}. Plot ({plot_residuals}).")
        else:
            RES = spspalinalg.tfqmr(A, b,
                                    x0=x0, tol=tol, maxiter=maxiter,
                                    atol=atol)

        x, info = RES
        assert x.__class__.__name__ == 'ndarray' and x.shape == (num_dofs,)

    else:
        info = None

    COMM.barrier()
    if RANK != MASTER_RANK:
        x = np.empty((num_dofs,), dtype=float)
    else:
        pass

    # noinspection PyUnboundLocalVariable
    COMM.Bcast([x, MPI.FLOAT], root=MASTER_RANK)
    info = COMM.bcast(info, root=MASTER_RANK)
    x = LocallyFullVector(x)
    t3 = time()
    message = f'scipy_sparse_linalg_tfqmr ({name}) = [SYSTEM-SHAPE: {A.shape}] >>> ' \
              f'costs {int((t3-t0)*100)/100}s with convergence info = {info}.'

    return x, info, 0, 0, message

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
