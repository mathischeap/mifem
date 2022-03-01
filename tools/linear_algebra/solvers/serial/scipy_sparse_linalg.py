# -*- coding: utf-8 -*-

"""
Here we store some serial_spspalinalg solvers.

"""
from root.config import *
from scipy import sparse as spspa
from scipy.sparse import linalg as spspalinalg
from tools.linear_algebra.data_structures.global_matrix.main import DistributedVector
from time import time


def spsolve(AA, bb):
    """

    :param AA:
    :param bb:
    :return: Return a tuple of 4 outputs:

            1. (DistributedVector) results -- The result vector.
            2. 0 --
            3. 0 --
            4. 0 --
            5. (str) message
    """
    cOmm.barrier()
    t0 = time()
    assert AA.__class__.__name__ == 'GlobalMatrix' and AA.IS_master_dominating, \
        " <serial> <scipy sparse linalg> <{spsolve}> A must be GlobalMatrix and master dominating."
    assert bb.__class__.__name__ == 'GlobalVector' and bb.IS_master_dominating, \
        " <serial> <scipy sparse linalg> <{spsolve}> b must be GlobalVector and master dominating."
    cOmm.barrier()
    t1 = time()
    # ...
    if rAnk == mAster_rank:
        shape = AA.shape[0]
        assert shape == bb.shape[0], f"A:{AA.shape} does not match b{bb.shape}."
        x = spspalinalg.spsolve(AA.M, bb.V)
        x = spspa.csr_matrix(x).T
    else:
        x = None
    cOmm.barrier()
    t2 = time()
    x = DistributedVector(x)
    assert x.IS_master_dominating
    cOmm.barrier()
    t3 = time()
    if rAnk == mAster_rank:
        message = f'scipy.sparse.linalg.spsolve {AA.shape}, ' \
                  f'cost {int((t2-t1)*100)/100}/{int((t3-t0)*100)/100}.'
    else:
        message = ''
    return x, 0, 0, 0, message


def gmres(*args, **kwargs):
    return ___scipy_sparse_linalg___("gmres", *args, **kwargs)

def gcrotmk(*args, **kwargs):
    return ___scipy_sparse_linalg___("gcrotmk", *args, **kwargs)

def bicgstab(*args, **kwargs):
    return ___scipy_sparse_linalg___("bicgstab", *args, **kwargs)


def ___scipy_sparse_linalg___(solver_name, AA, bb, X0, restart=100, maxiter=1000, tol=1e-5):
    """

    :param AA:
    :param bb:
    :param X0:
    :param restart:
    :param maxiter:
    :param tol:
    :return: Return a tuple of 4 outputs:

            1. (DistributedVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * <0 : illegal input or breakdown

            3. (float) beta -- The residual.
            4. (int) ITER -- The number of outer iterations.
            5. (str) message
    """
    cOmm.barrier()
    t0 = time()
    assert AA.__class__.__name__ == 'GlobalMatrix', f" <serial> <scipy sparse linalg> <{solver_name}> A must be GlobalMatrix"
    assert bb.__class__.__name__ == 'GlobalVector', f" <serial> <scipy sparse linalg> <{solver_name}> b must be GlobalVector"
    assert X0.__class__.__name__ == 'DistributedVector', f" <serial> <scipy sparse linalg> <{solver_name}> x0 must be DistributedVector"
    assert maxiter >= 1, f" <serial> <scipy sparse linalg> <{solver_name}> maxiter must be >= 1."
    assert restart >= 3, f" <serial> <scipy sparse linalg> <{solver_name}> restart must be >= 3."
    assert tol > 0, f" <serial> <scipy sparse linalg> <{solver_name}> tol must be > 0."


    assert AA.__class__.__name__ == 'GlobalMatrix' and AA.IS_master_dominating, \
        f" <serial> <scipy sparse linalg> <{solver_name}> A must be GM and master dominating."
    assert bb.__class__.__name__ == 'GlobalVector' and bb.IS_master_dominating, \
        f" <serial> <scipy sparse linalg> <{solver_name}> b must be GM and master dominating."

    if X0.IS_master_dominating:
        if rAnk == mAster_rank:
            x0 = X0.V.T.toarray()[0]
    else:
        x0 = X0.___PRIVATE_gather_V_to_core___(clean_local=True)

    cOmm.barrier()
    t1 = time()
    # ...
    if rAnk == mAster_rank:
        # preconditioning
        invM = Jacobian_preconditioner(AA.M)


        A = AA.M
        if bb.V.__class__.__name__ == 'ndarray':
            b = bb.V
        else:
            b = bb.V.toarray().ravel()


        solver = getattr(spspalinalg, solver_name)
        if solver_name == 'gmres':
            x, info = solver(A, b, x0=x0, M=invM, maxiter=maxiter, tol=tol, atol=tol,
                restart=restart)
        elif solver_name in ('bicgstab', 'gcrotmk'):
            x, info = solver(AA.M, b, x0=x0, M=invM, maxiter=maxiter, tol=tol, atol=tol)
        else:
            raise NotImplementedError(f"no solver named {solver_name} implemented.")

        x = spspa.csr_matrix(x).T
    else:
        x = None
        info = None
    cOmm.barrier()
    t2 = time()

    info = cOmm.bcast(info, root=mAster_rank)
    x = DistributedVector(x)
    assert x.IS_master_dominating

    cOmm.barrier()
    t3 = time()

    if rAnk == mAster_rank:
        message = f'scipy.sparse.linalg.{solver_name} {AA.shape}, ' \
                  f'cost {int((t2-t1)*100)/100}/{int((t3-t0)*100)/100}, ' \
                  f'[convergence_info: {info}], tol={tol}, maxiter={maxiter}.'
    else:
        message = ''
    if solver_name == 'gmres': message = message[:-1] + f', restart={restart}.'

    return x, info, tol, info, message



def Jacobian_preconditioner(A):
    """Produce a Jacobian_preconditioner for sparse matrix A.

    :param A:
    :return:
    """
    diag = A.diagonal()
    diag[diag==0] = 1.0
    diag[0] = float(diag[0]) # make sure reciprocal works correctly.
    diag = np.reciprocal(diag)
    invA = spspa.dia_matrix((diag, [0,]), shape=A.shape)
    return invA
