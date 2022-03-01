"""
BiCGSTAB

Not working.

"""



from root.config import *
from tools.linear_algebra.preconditioners.allocator import PreconditionerAllocator
from screws.miscellaneous import MyTimer
from tools.linear_algebra.data_structures.global_matrix.main import LocallyFullVector
from screws.exceptions import LinerSystemSolverDivergenceError



from tools.linear_algebra.solvers.parallel.GMRES import ___gmres_stop_criterion___



def solve(A, b, x0, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=(None, dict()), COD=True, routine='auto'):
    """

    :param A: GlobalMatrix
    :param b: GlobalVector
    :param x0: LocallyFullVector
    :param maxiter:
    :param tol: tolerance.
    :param atol: absolute tolerance.
    :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
    :param COD: Clear Original Data?
    :param routine: Which particular routine we are going to use?
    :return: Return a tuple of 5 outputs:

            1. (DistributedVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * -1 : divergence


            3. (float) beta -- The residual.
            4. (int) ITER -- The number of outer iterations.
            5. (str) message

    """


    message = "BiCGSTAB-" + MyTimer.current_time()
    assert A.__class__.__name__ == "GlobalMatrix", \
                     f"A needs to be a GlobalMatrix'. Now I get {A.__class__}."

    assert b.__class__.__name__ == "GlobalVector", \
                     f"b needs to be a 'GlobalVector'. Now I get {b.__class__}."

    assert x0.__class__.__name__ == "LocallyFullVector", \
                     f"x0 needs to be a 'LocallyFullVector'. Now I get {b.__class__}."
    assert maxiter >= 1 and maxiter % 1 == 0, f"maxiter={maxiter} must be >= 1."
    assert tol > 0 and atol > 0, f"tol={tol} and atol={atol} wrong, they must be > 0."

    # -------  Decide preconditioner -----------------------------------------------------------------------------------
    if preconditioner is None: preconditioner = (None, dict())

    preconditioner_ID, preconditioner_kwargs = preconditioner
    if preconditioner_ID is not None:
        assert preconditioner_ID in PreconditionerAllocator.___defined_preconditioners___(), \
            f"preconditioner={preconditioner_ID} is not coded, try one of " \
            f"{PreconditionerAllocator.___defined_preconditioners___().keys()}"
        preconditioner = PreconditionerAllocator(preconditioner_ID)(A, **preconditioner_kwargs)
    else:
        preconditioner = None

    # -------  Decide routine ------------------------------------------------------------------------------------------
    if routine == 'auto':
        ROUTINE = ___mpi_v0_BiCGSTAB___ # in the future, we may want to make de function to decide which one is the best for particular matrices.
    else:
        if routine == 'mpi_v0':
            ROUTINE = ___mpi_v0_BiCGSTAB___
        else:
            raise Exception(f"routine={routine} is wrong.")

    # ---------- Do the computation ------------------------------------------------------------------------------------
    results, info, beta, ITER, solver_message = \
    ROUTINE(A, b, x0, maxiter=maxiter, tol=tol, atol=atol, preconditioner=preconditioner, COD=COD)

    MESSAGE =  message + '-' + solver_message
    #===================================================================================================================

    return results, info, beta, ITER, MESSAGE


def ___mpi_v0_BiCGSTAB___(lhs, rhs, X0, maxiter=3, tol=1e-3, atol=1e-4, preconditioner=None, COD=True):
    """

    :param lhs: GlobalMatrix
    :param rhs: LocallyFullVector
    :param X0: LocallyFullVector
    :param maxiter:
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

        if applying_method == 'left_multiply':
            invM = preconditioner.M
            A = invM @ A
            f = invM @ f
        else:
            raise NotImplementedError(f"We did not yet code preconditioning for the "
                                      f"routine: ___mpi_v2_gmres___ using method: <{applying_method}>.")


    if COD:  # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_ = None

    if rAnk == mAster_rank:
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
        cOmm.Reduce(_r0_acd_, r0, op=MPI.SUM, root=mAster_rank)

        if rAnk == mAster_rank:
            if ITER == 0: hr0 = r0
            beta = np.sum(r0 ** 2) ** 0.5
        else:
            beta = 0

        beta = cOmm.bcast(beta, root=mAster_rank)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,] # this is right, do not initial BETA as an empty list.
        if len(BETA) > 20: BETA = BETA[:1] + BETA[-2:]
        BETA.append(beta)
        JUDGE, stop_iteration, info, JUDGE_explanation = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        ITER += 1

        if rAnk == mAster_rank and ITER > 1:
            r0 = s - omega1 * t
            rho0 = rho1
            omega0 = omega1
            p0 = p1
            v0 = v1

        if rAnk == mAster_rank:
            rho1 = np.sum(hr0 * r0)
            beta = (rho1/rho0) * (alpha/omega0)
            p1 = r0 + beta * (p0 - omega0*v0)

        cOmm.Bcast([p1, MPI.FLOAT], root=mAster_rank)

        _v1_acd_ = A @ p1

        cOmm.Reduce(_v1_acd_, v1, op=MPI.SUM, root=mAster_rank)

        if rAnk == mAster_rank:
            alpha = rho1 / np.sum(hr0 * v1)
            s = r0 - alpha * v1

        cOmm.Bcast([s, MPI.FLOAT], root=mAster_rank)

        _t_acd_ = A @ s
        cOmm.Reduce(_t_acd_, t, op=MPI.SUM, root=mAster_rank)

        if rAnk == mAster_rank:
            omega1 = np.sum(t * s) / np.sum(t ** 2)
            x0 += alpha * p1 + omega1 * s


        cOmm.Bcast([x0, MPI.FLOAT], root=mAster_rank)


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