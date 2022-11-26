# -*- coding: utf-8 -*-

from root.config.main import *
from tools.linearAlgebra.dataStructures.globalMatrix.main import LocallyFullVector
from components.exceptions import LinerSystemSolverDivergenceError
from tools.linearAlgebra.solvers.regular.GMRES.helpers.components.stop_criterion import ___gmres_stop_criterion___


def ___mpi_v0_BiCGSTAB___(lhs, rhs, X0, maxiter=3, tol=1e-3, atol=1e-4, preconditioner=None, COD=True):
    """

    :param lhs: GlobalMatrix
    :param rhs: LocallyFullVector
    :param X0: LocallyFullVector
    :param maxiter: int, str
        A positive integer.

        if maxiter is a str, it must be a numeric str, and it means it is a
        strong maxiter, that is no matter what happened, we will iterate the
        solver for this many times. So it is a forced amount of iterations.

    :param tol: relative tolerance.
    :param atol: absolute tolerance.
    :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
    :param COD: Clear Original Data?
    :return: Return a tuple of 5 outputs:

            1. (LocallyFullVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : (convergence to tolerance not achieved) `Info` means number of iterations
                * <0 : illegal input or breakdown

            3. (float) beta -- The residual.
            4. (int) ITER -- The number of outer iterations.
            5. (str) message

    """

    Time_start = MPI.Wtime()

    A = lhs.M
    f = rhs.V.toarray().ravel()
    x0 = np.array(X0.V) # just to make another copy of x0: IMPORTANT TO do SO!

    shape0, shape1 = A.shape
    assert f.shape[0] == x0.shape[0] == shape0 == shape1, "Ax=f shape dis-match."

    if preconditioner is not None:
        applying_method = preconditioner.applying_method

        if applying_method == 'left_multiply_invM':
            invM = preconditioner.invM
            A = invM @ A
            f = invM @ f
        else:
            raise NotImplementedError(f"We did not yet code preconditioning for the "
                                      f"routine: ___mpi_v1_gmres___ using method: <{applying_method}>.")


    if COD:  # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_ = None

    if RANK == MASTER_RANK:
        r0 = np.empty((shape0,), dtype=float)
        rho0   = 1
        omega0 = 1
        p0 = np.zeros((shape0,), dtype=float)
        v0 = np.zeros((shape0,), dtype=float)

        alpha = 1
        v1 = np.empty((shape0,), dtype=float)
        t  = np.empty((shape0,), dtype=float)

    else:
        r0 = None
        t = None
        v1 = None

        p1 = np.zeros((shape0,), dtype=float)
        s  = np.empty((shape0,), dtype=float)

    ITER = 0
    BETA = None

    while 1:
        _r0_acd_ = f - A @ x0 # acd stands for `all cores distributed`
        COMM.Reduce(_r0_acd_, r0, op=MPI.SUM, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            if ITER == 0: hr0 = r0
            beta = np.sum(r0 ** 2) ** 0.5
        else:
            beta = 0

        beta = COMM.bcast(beta, root=MASTER_RANK)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,] # this is right, do not initial BETA as an empty list.
        if len(BETA) > 20: BETA = BETA[:1] + BETA[-2:]
        BETA.append(beta)
        JUDGE, stop_iteration, info, JUDGE_explanation = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        ITER += 1

        if RANK == MASTER_RANK and ITER > 1:
            r0 = s - omega1 * t
            rho0 = rho1
            omega0 = omega1
            p0 = p1
            v0 = v1

        if RANK == MASTER_RANK:
            rho1 = np.sum(hr0 * r0)
            beta = (rho1/rho0) * (alpha/omega0)
            p1 = r0 + beta * (p0 - omega0*v0)

        COMM.Bcast([p1, MPI.FLOAT], root=MASTER_RANK)

        _v1_acd_ = A @ p1

        COMM.Reduce(_v1_acd_, v1, op=MPI.SUM, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            alpha = rho1 / np.sum(hr0 * v1)
            s = r0 - alpha * v1

        COMM.Bcast([s, MPI.FLOAT], root=MASTER_RANK)

        _t_acd_ = A @ s
        COMM.Reduce(_t_acd_, t, op=MPI.SUM, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            omega1 = np.sum(t * s) / np.sum(t ** 2)
            x0 += alpha * p1 + omega1 * s


        COMM.Bcast([x0, MPI.FLOAT], root=MASTER_RANK)


    x0 = LocallyFullVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    Time_end = MPI.Wtime()

    COST_total = Time_end - Time_start
    message = f" mpi_v0_BiCGSTAB = [SYSTEM-SHAPE: {A.shape}] [ITER={ITER}][residual=%.2e] costs %.2f, " \
              f"convergence info={info}, maxiter={maxiter}, " \
              f"stop_judge={JUDGE}: {JUDGE_explanation}]"%(beta, COST_total)

    return x0, info, beta, ITER, message