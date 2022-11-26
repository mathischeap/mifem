# -*- coding: utf-8 -*-
from root.config.main import *
from tools.linearAlgebra.dataStructures.globalMatrix.main import LocallyFullVector
from components.exceptions import LinerSystemSolverDivergenceError
from tools.linearAlgebra.solvers.regular.GMRES.helpers.components.stop_criterion import ___gmres_stop_criterion___
from tools.linearAlgebra.solvers.regular.GMRES.helpers.components.residual_ploter import ___gmres_plot_residuals___

def ___mpi_v0_gmres___(lhs, rhs, X0, restart=100, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=None,
                       COD=True, name=None, plot_residuals=False):
    """
    :param lhs: GlobalMatrix
    :param rhs: GlobalVector
    :param X0: LocallyFullVector
    :param restart:
    :param maxiter: int, str
        A positive integer.

        if maxiter is a str, it must be a numeric str, and it means it is a
        strong maxiter, that is no matter what happened, we will iterate the
        solver for this many times. So it is a forced amount of iterations.

    :param tol: relative tolerance.
    :param atol: absolute tolerance.
    :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
    :param COD: Clear Original Data?
    :param name: The name of this solving process.
    :param plot_residuals: bool
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
                                      f"routine: ___mpi_v0_gmres___ using method: <{applying_method}>.")

    if COD: # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_ = None

    ITER = 0
    BETA = None

    if RANK == MASTER_RANK:
        AVJ = np.empty((shape0,), dtype=float)
        Hm = np.zeros((restart + 1, restart), dtype=float)
        VV = np.empty((shape0,), dtype=float)
        Vm = np.empty((restart, shape0), dtype=float)

        if plot_residuals:
            residuals = list()

    else:
        AVJ = None
        VV = None

    while 1: # always do till break.

        v0 = f - A @ x0

        COMM.Reduce(v0, VV, op=MPI.SUM, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            beta = np.sum(VV ** 2) ** 0.5
            v0 = VV / beta
        else:
            beta = None

        COMM.Bcast([v0, MPI.FLOAT], root=MASTER_RANK)
        beta = COMM.bcast(beta, root=MASTER_RANK)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,] # this is right, do not initial BETA as an empty list.
        if len(BETA) > 20: BETA = BETA[:1] + BETA[-2:]
        BETA.append(beta)
        if RANK == MASTER_RANK:
            if plot_residuals:
                # noinspection PyUnboundLocalVariable
                residuals.append(beta)
        JUDGE, stop_iteration, info, JUDGE_explanation = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        if RANK == MASTER_RANK:
            # noinspection PyUnboundLocalVariable
            Vm[0] = v0
        else:
            Vm = v0

        for j in range(restart):
            if RANK == MASTER_RANK:
                Avj = A @ Vm[j]
            else:
                Avj = A @ Vm

            COMM.Reduce(Avj, AVJ, op=MPI.SUM, root=MASTER_RANK)

            if RANK == MASTER_RANK:

                sum_Hij_vi = 0

                for i in range(j+1):
                    Hij = np.sum(AVJ * Vm[i])
                    # noinspection PyUnboundLocalVariable
                    Hm[i,j] = Hij
                    sum_Hij_vi += Hij * Vm[i]

                hat_v_jp1 = AVJ - sum_Hij_vi
                Hm[j + 1, j] = np.sum(hat_v_jp1**2) ** 0.5

                if j < restart - 1:
                    v_jp1 = hat_v_jp1 / Hm[j+1, j]
                else:
                    del v_jp1

            else:
                if j == 0:
                    v_jp1 = np.empty((shape0,), dtype=float)
                else:
                    pass

            if j < restart - 1:
                # noinspection PyUnboundLocalVariable
                COMM.Bcast([v_jp1, MPI.FLOAT], root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    Vm[j+1] = v_jp1
                else:
                    Vm = v_jp1

        if RANK == MASTER_RANK:
            HmT = Hm.T
            ls_A = HmT @ Hm
            ls_b = HmT[:,0] * beta
            ym = np.linalg.solve(ls_A, ls_b)
            del HmT, ls_A, ls_b
            x0 += ym.T @ Vm

        COMM.Bcast([x0, MPI.FLOAT], root=MASTER_RANK)
        ITER += 1

    x0 = LocallyFullVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    Time_end = MPI.Wtime()

    COST_total = Time_end - Time_start
    message = f" mpi_v0_gmres = [SYSTEM-SHAPE: {A.shape}] [ITER={ITER}] " \
              f"[residual=%.2e] costs %.2f, " \
              f"convergence info={info}, restart={restart}, maxiter={maxiter}, " \
              f"stop_judge={JUDGE}: {JUDGE_explanation}]"%(beta, COST_total)

    if RANK == MASTER_RANK:
        if plot_residuals:
            ___gmres_plot_residuals___(np.array(residuals), name, 'mpi_gmres_v1')

    return x0, info, beta, ITER, message