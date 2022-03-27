




from root.config.main import *
from tools.linear_algebra.data_structures.global_matrix.main import LocallyFullVector
from screws.exceptions import LinerSystemSolverDivergenceError

from tools.linear_algebra.solvers.parallel.GMRES.helpers.components.stop_criterion import ___gmres_stop_criterion___
from tools.linear_algebra.solvers.parallel.GMRES.helpers.components.residual_ploter import ___gmres_plot_residuals___




def ___mpi_v0_LGMRES___(lhs, rhs, X0, m=100, k=10, maxiter=50, tol=1e-3, atol=1e-4, preconditioner=None,
                        COD=True, name=None, plot_residuals=False):
    """

    :param lhs: GlobalMatrix
    :param rhs: GlobalVector
    :param X0: LocallyFullVector
    :param m:
    :param k:
    :param maxiter:
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
    _m_, _k_ = m, k
    restart = m + k

    local_ind = list()
    for i in range(restart):
        if (i % sIze) == rAnk:
            local_ind.append(i)

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
                                      f"routine: ___mpi_v1_gmres___ using method: <{applying_method}>.")

    if COD: # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_  = None

    ITER = 0
    BETA = None

    AVJ = np.empty((shape0,), dtype=float)
    Hm = np.zeros((restart + 1, restart), dtype=float) # In the future, we can change this to sparse matrix.
    if rAnk == mAster_rank:
        VV = np.empty((shape0,), dtype=float)
        Vm = np.empty((restart, shape0), dtype=float)

        SUM_Hij_vi = np.empty((shape0,), dtype=float)
        Vs = None
        local_ind_dict = None

        if plot_residuals:
            residuals = list()
    else:
        local_ind_dict = dict()
        for i, ind in enumerate(local_ind):
            local_ind_dict[ind] = i

        VV = None
        SUM_Hij_vi = None
        Vs = np.empty((len(local_ind), shape0), dtype=float)
        Vm = None

    if rAnk == sEcretary_rank:
        HM = np.empty((restart + 1, restart), dtype=float)  # use to combine all Hm.
    else:
        HM = None

    ZZZ = dict()
    AZ_cache = dict()

    while 1: # always do till break.

        if ITER < _k_:
            k = ITER
            m = restart - k
        else:
            m = _m_
            k = _k_

        v0 = f - A @ x0
        cOmm.Reduce(v0, VV, op=MPI.SUM, root=mAster_rank)
        if rAnk == mAster_rank:
            beta = np.sum(VV ** 2) ** 0.5
            v0 = VV / beta
        else:
            beta = None

        cOmm.Bcast([v0, MPI.FLOAT], root=mAster_rank)
        beta = cOmm.bcast(beta, root=mAster_rank)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,] # this is right, do not initial BETA as an empty list.
        if len(BETA) > 20: BETA = BETA[:1] + BETA[-2:]
        BETA.append(beta)
        if rAnk == mAster_rank:
            if plot_residuals:
                # noinspection PyUnboundLocalVariable
                residuals.append(beta)

        JUDGE, stop_iteration, info, JUDGE_explanation = \
            ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        if rAnk == mAster_rank:
            Vm[0] = v0
        else:
            Vm = v0
            if 0 in local_ind:
                Vs[local_ind_dict[0]] = v0

        for j in range(restart):

            if j < m:
                if rAnk == mAster_rank:
                    Avj = A @ Vm[j]
                else:
                    Avj = A @ Vm
                    if j in local_ind: Vs[local_ind_dict[j]] = Vm
            else:
                index = ITER + m - 1 - j

                if index in AZ_cache:
                    pass
                else:
                    if ITER > _k_: del AZ_cache[ITER - _k_ - 1]
                    AZ_cache[index] = A @ ZZZ[index]

                Avj = AZ_cache[index]

                if rAnk == mAster_rank:
                    pass
                else:
                    if j in local_ind: Vs[local_ind_dict[j]] = Vm

            cOmm.Allreduce(Avj, AVJ, op=MPI.SUM)

            sum_Hij_vi = np.zeros((shape0,), dtype=float)

            if rAnk == mAster_rank:
                for i in local_ind:
                    if i <= j:
                        _ = Vm[i]
                        Hij = np.sum(AVJ * _)
                        sum_Hij_vi += Hij * _
                        Hm[i,j] = Hij
                    else:
                        break
            else:
                for i in local_ind:
                    if i <= j:
                        _ = Vs[local_ind_dict[i]]
                        Hij = np.sum(AVJ * _)
                        sum_Hij_vi += Hij * _
                        Hm[i,j] = Hij
                    else:
                        break

            cOmm.Reduce(sum_Hij_vi, SUM_Hij_vi, op=MPI.SUM, root=mAster_rank)

            if rAnk == mAster_rank:
                hat_v_jp1 = AVJ - SUM_Hij_vi
                Hm[j + 1, j] = np.sum(hat_v_jp1 ** 2) ** 0.5

                if j < restart - 1:
                    v_jp1 = hat_v_jp1 / Hm[j + 1, j]
                else:
                    del v_jp1

            else:
                if j == 0:
                    v_jp1 = np.empty((shape0,), dtype=float)
                else:
                    pass

            if j < restart - 1:
                # noinspection PyUnboundLocalVariable
                cOmm.Bcast([v_jp1, MPI.FLOAT], root=mAster_rank)
                if rAnk == mAster_rank:
                    Vm[j + 1] = v_jp1
                else:
                    Vm = v_jp1

        cOmm.Reduce(Hm, HM, op=MPI.SUM, root=sEcretary_rank)

        if rAnk == sEcretary_rank:
            HMT = HM.T
            ls_A = HMT @ HM
            ls_b = HMT[:,0] * beta
            ym = np.linalg.solve(ls_A, ls_b)
        else:
            pass

        if rAnk == mAster_rank:
            if k == 0:
                Ws = Vm
            else:
                DL = list()
                iL = [ITER - _ - 1 for _ in range(k)]
                for _ in iL: DL.append(ZZZ[_])
                Ws = np.vstack((Vm[:m], np.array(DL)))
        else:
            pass

        if mAster_rank != sEcretary_rank:
            if rAnk == mAster_rank:
                ym = np.empty((restart,), dtype=float)
                cOmm.Recv([ym, MPI.FLOAT], source=sEcretary_rank, tag=ITER)
            elif rAnk == sEcretary_rank:
                # noinspection PyUnboundLocalVariable
                cOmm.Send([ym, MPI.FLOAT], dest=mAster_rank, tag=ITER)
            else:
                pass
        else:
            pass

        if rAnk == mAster_rank:
            # noinspection PyUnboundLocalVariable
            ym = ym.T @ Ws
        else:
            ym = np.empty((shape0,), dtype=float) # Important: renew ``ym`` every single iteration.

        cOmm.Bcast([ym, MPI.FLOAT], root=mAster_rank)

        if ITER >= _k_ > 0: del ZZZ[ITER-_k_]
        ZZZ[ITER] = ym

        x0 += ym
        ITER += 1

    x0 = LocallyFullVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    Time_end = MPI.Wtime()

    COST_total = Time_end - Time_start
    message = f" mpi_v0_LGMRES = [SYSTEM-SHAPE: {A.shape}] [ITER={ITER}][residual=%.2e] costs %.2f, " \
              f"convergence info={info}, m={_m_}, k={_k_} maxiter={maxiter}, " \
              f"stop_judge={JUDGE}: {JUDGE_explanation}]"%(beta, COST_total)



    if rAnk == mAster_rank:
        if plot_residuals:
            ___gmres_plot_residuals___(np.array(residuals), name, 'mpi_LGMRES_v0')

    return x0, info, beta, ITER, message