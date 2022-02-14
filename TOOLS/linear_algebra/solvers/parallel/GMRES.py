"""
See paper:

[Saad and Schultz, GMRES: A GENERALIZED MINIMAL RESIDUAL ALGORITHM FOR SOLVING NON-SYMMETRIC LINEAR SYSTEMS,
SIAM J. ScI. STAT. COMPUT., 1986]

"""

import os

from root.config import *
from TOOLS.linear_algebra.preconditioners.allocator import PreconditionerAllocator
from TOOLS.linear_algebra.data_structures import LocallyFullVector
from SCREWS.exceptions import LinerSystemSolverDivergenceError
from SCREWS.miscellaneous import MyTimer

import matplotlib.pyplot as plt

def solve(A, b, x0, restart=100, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=(None, dict()), COD=True,
    routine='auto', loading_factor = 3e6, name=None, plot_residuals=False, **kwargs):
    """

    :param A: GlobalMatrix
    :param b: GlobalVector
    :param x0: LocallyFullVector
    :param restart:
    :param maxiter:
    :param tol: tolerance.
    :param atol: absolute tolerance.
    :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
    :param COD: Clear Original Data?
    :param routine: Which particular routine we are going to use?
    :param loading_factor: A factor that decide the limit a single core (the master core) is going to handle. When the
        computational power of a single core (the master core) increases, we can make this factor larger.
    :param name: The name of this solving process.
    :param plot_residuals: bool, if we plot the residuals.
    :param kwargs: possible other args for particular routine.
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

    message = "GMRES-" + MyTimer.current_time()
    assert str(A.__class__) == "<class 'TOOLS.linear_algebra.data_structures.GlobalMatrix'>", \
                     f"A needs to be a 'TOOLS.linear_algebra.data_structures.GlobalMatrix'. Now I get {A.__class__}."

    assert str(b.__class__) == "<class 'TOOLS.linear_algebra.data_structures.GlobalVector'>", \
                     f"b needs to be a 'TOOLS.linear_algebra.data_structures.GlobalVector'. Now I get {b.__class__}."

    assert str(x0.__class__) == "<class 'TOOLS.linear_algebra.data_structures.LocallyFullVector'>", \
                     f"x0 needs to be a 'TOOLS.linear_algebra.data_structures.LocallyFullVector'. Now I get {b.__class__}."

    assert maxiter >= 1 and maxiter % 1 == 0, f"maxiter={maxiter} must be >= 1."
    assert restart >= 3 and restart % 1 == 0, f"restart={restart} must be >= 3."
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
        shape0, _ = A.shape

        if shape0 * restart < loading_factor: # when the total loading is less than the loading_factor,
                                              # we use routine v0 which is faster when the loading is low because it
                                              # needs less communications
            ROUTINE = ___mpi_v0_gmres___
        else: # otherwise, we use routine v2 which is faster when the loading is high.
            ROUTINE = ___mpi_v2_gmres___

    else:
        if routine == 'mpi_v0':
            ROUTINE = ___mpi_v0_gmres___
        # elif routine == 'mpi_v1': # routine v1 does not work...
        #     ROUTINE = ___mpi_v1_gmres___
        elif routine == 'mpi_v2':
            ROUTINE = ___mpi_v2_gmres___
        else:
            raise Exception(f"routine={routine} is wrong.")

    # ---------- Do the computation ------------------------------------------------------------------------------------
    results, info, beta, ITER, solver_message = \
    ROUTINE(A, b, x0, restart=restart, maxiter=maxiter, tol=tol, atol=atol, preconditioner=preconditioner,
            COD=COD, name=name, plot_residuals=plot_residuals)

    _ = kwargs # trivial; just leave freedom for future updates for kwargs.

    MESSAGE =  message + '-' + solver_message
    #===================================================================================================================

    return results, info, beta, ITER, MESSAGE



def ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA):
    """
    :param tol: relative tolerance.
    :param atol: absolute tolerance
    :param ITER:
    :param maxiter:
    :param BETA: A list of beta (residual) of some recent iterations.
    :return:
    """
    assert tol < 0.01, f"tol={tol} too large, should be < 0.01."

    # noinspection PyUnusedLocal
    info = 'TBD'

    beta0 = BETA[0]
    beta = BETA[-1]
    judge_1 = beta < atol # judge 1: reach absolute tolerance.
    judge_2 = ITER >= maxiter # judge 2: reach max iteration number
    # judge 3: divergence
    if BETA[-1] > BETA[-2]: # error grows after one iteration
        if BETA[-2] > 1 and (BETA[-1]-BETA[-2]) > 100 * BETA[-2]:
            judge_3 = True
        elif BETA[-1] > 10e6:
            judge_3 = True
        elif (BETA[-1]-BETA[-2]) > 100:
            judge_3 = True
        else:
            judge_3 = False
    else:
        judge_3 = False

    # judge 4: reach relative tol.
    if beta < beta0:
        progress = beta0 - beta
        if progress / beta0 < tol: # reach relative tol.
            judge_4 = True
        else:
            judge_4 = False
    else:
        judge_4 = False

    # judge_5: slow converging
    beta_old = BETA[-2]
    if beta < beta_old:
        progress = beta_old - beta
        if progress / beta_old < tol: # slow converging
            judge_5 = True
        else:
            judge_5 = False
    else:
        judge_5 = False

    # ...

    if judge_1 or judge_2 or judge_3 or judge_4 or judge_5:

        stop_iteration = True


        if judge_1: # reach atol
            info = 0
            JUDGE = 1
            JUDGE_explanation = 'reach absolute tol'

        elif judge_2: # reach maxiter
            info = ITER
            JUDGE = 2
            JUDGE_explanation = 'reach maxiter'


        elif judge_3: # diverging
            info = -1
            JUDGE = 3
            JUDGE_explanation = 'diverging'


        elif judge_4: # reach tol
            info = 0
            JUDGE = 4
            JUDGE_explanation = 'reach relative tol'


        elif judge_5: # very slow converging; the progress is lower than the tol
            info = ITER
            JUDGE = 5
            JUDGE_explanation = 'very slow converging'


        else:
            raise Exception()


    else: # do not stop iterations.
        stop_iteration = False
        info = None
        JUDGE = 0
        JUDGE_explanation = ''

    assert stop_iteration in (True, False), "stop_iteration has to be set."
    assert info != 'TBD', "info has to be updated"

    return JUDGE, stop_iteration, info, JUDGE_explanation




def ___gmres_plot_residuals___(residuals, solve_name, scheme_name):
    """

    :param residuals: We plot these residuals.
    :param solve_name: We will save the plot using this name.
    :return:
    """
    assert solve_name is not None, f"to plot residuals, give the solving process a name."

    fig = plt.figure()
    plt.semilogy(residuals)
    plt.title(scheme_name + ': ' + solve_name  + '@' + MyTimer.current_time())
    plt.ylabel('residuals')
    plt.xlabel('iterations')

    if os.path.isdir("__RESIDUALS__"):
        pass
    else:
        os.mkdir('__RESIDUALS__')

    plt.savefig(f'./__RESIDUALS__/'
        f'{scheme_name + "_" + solve_name + "_" + MyTimer.current_time_with_no_special_characters()}.png')

    plt.close(fig)




def ___mpi_v0_gmres___(lhs, rhs, X0, restart=100, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=None,
                       COD=True, name=None, plot_residuals=False):
    """

    :param lhs: GlobalMatrix
    :param rhs: GlobalVector
    :param X0: LocallyFullVector
    :param restart:
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

    Time_start = MPI.Wtime()
    A = lhs.M
    f = rhs.V.toarray().ravel()
    x0 = np.array(X0.V) # just to make another copy of x0: IMPORTANT TO DO SO!

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
                                      f"routine: ___mpi_v0_gmres___ using method: <{applying_method}>.")


    if COD: # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_ = None

    ITER = 0
    BETA = None

    if rAnk == mAster_rank:
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
                residuals.append(beta)
        JUDGE, stop_iteration, info, JUDGE_explanation = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        if rAnk == mAster_rank:
            # noinspection PyUnboundLocalVariable
            Vm[0] = v0
        else:
            Vm = v0

        for j in range(restart):
            if rAnk == mAster_rank:
                Avj = A @ Vm[j]
            else:
                Avj = A @ Vm

            cOmm.Reduce(Avj, AVJ, op=MPI.SUM, root=mAster_rank)

            if rAnk == mAster_rank:

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
                cOmm.Bcast([v_jp1, MPI.FLOAT], root=mAster_rank)
                if rAnk == mAster_rank:
                    Vm[j+1] = v_jp1
                else:
                    Vm = v_jp1

        if rAnk == mAster_rank:
            HmT = Hm.T
            ls_A = HmT @ Hm
            ls_b = HmT[:,0] * beta
            ym = np.linalg.solve(ls_A, ls_b)
            del HmT, ls_A, ls_b
            x0 += ym.T @ Vm

        cOmm.Bcast([x0, MPI.FLOAT], root=mAster_rank)
        ITER += 1

    x0 = LocallyFullVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    Time_end = MPI.Wtime()

    COST_total = Time_end - Time_start
    message = f" mpi_v0_gmres = [SYSTEM-SHAPE: {A.shape}] [ITER={ITER}][residual=%.2e] costs %.2f, " \
              f"convergence info={info}, restart={restart}, maxiter={maxiter}, " \
              f"stop_judge={JUDGE}: {JUDGE_explanation}]"%(beta, COST_total)

    if rAnk == mAster_rank:
        if plot_residuals:
            ___gmres_plot_residuals___(np.array(residuals), name, 'mpi_gmres_v1')

    return x0, info, beta, ITER, message



# def ___mpi_v1_distributor___(loading_factor, shape0, restart):
#     for i in range(2, sIze + 1):
#         loading_number = loading_factor * 1e8
#         total_number = shape0 * restart
#         if total_number < loading_number:
#
#             Vm_Cor = [mAster_rank, ]
#             Vm_Dis = [restart, ]
#             Vm_Ind = [range(0, restart), ]
#
#         elif (total_number / loading_number) >= sIze:  # then all cores will be used
#
#             if shape0 <= (loading_number / 5):
#                 num = restart
#                 parts = sIze
#                 distribution = [num // parts + (1 if x < num % parts else 0) for x in range(parts)]
#                 Vm_Dis = distribution
#                 other_cores = list()
#                 for i in range(sIze):
#                     if i != mAster_rank and i != sEcretary_rank:
#                         other_cores.append(i)
#                 Vm_Cor = [mAster_rank, sEcretary_rank] + other_cores
#                 distribution = [0, ] + distribution
#                 Vm_Ind = [range(sum(distribution[:i + 1]), sum(distribution[:i + 2])) for i in range(sIze)]
#
#             else:
#                 raise NotImplementedError("For large #DOF's.")
#
#         else:
#
#             Num_Cores = int(total_number / loading_number)
#             assert Num_Cores >= 1 and Num_Cores< sIze, "must be!"
#             if Num_Cores == 1:
#                 Num_Cores = 2
#             num = restart
#             parts = Num_Cores
#             distribution = [num // parts + (1 if x < num % parts else 0) for x in range(parts)]
#             Vm_Dis = distribution
#             if parts == 2:
#                 Vm_Cor = [mAster_rank, sEcretary_rank]
#             else:
#                 other_cores = list()
#                 for i in range(sIze):
#                     if i != mAster_rank and i != sEcretary_rank:
#                         other_cores.append(i)
#
#                 Vm_Cor = [mAster_rank, sEcretary_rank] + other_cores[:Num_Cores-2]
#
#             distribution = [0, ] + distribution
#             Vm_Ind = [range(sum(distribution[:i + 1]), sum(distribution[:i + 2])) for i in range(sIze)]
#
#
#
#         if saFe_mode:  # do some regular checks.
#             assert len(set(Vm_Cor)) == len(Vm_Cor) and mAster_rank == Vm_Cor[0]
#             if len(Vm_Cor) > 1 and sEcretary_rank != mAster_rank: assert sEcretary_rank == Vm_Cor[1]
#             assert Vm_Ind[-1].stop == restart
#             assert Vm_Ind[0].start == 0
#             assert sum(Vm_Dis) == restart and len(Vm_Dis) == len(Vm_Cor) == len(Vm_Ind)
#             for i, vmi in enumerate(Vm_Ind[1:-1]):
#                 vmi_m1, vmi_p1 = Vm_Ind[i], Vm_Ind[i + 2]
#                 assert vmi.start == vmi_m1.stop
#                 assert vmi.stop == vmi_p1.start
#
#         if rAnk == mAster_rank: print(Vm_Cor, Vm_Dis, Vm_Ind)
#
#         return Vm_Cor, Vm_Dis, Vm_Ind
#
# def ___mpi_v1_gmres___(lhs, rhs, X0, restart=100, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=None, COD=True, loading_factor=1):
#     """
#     In this version, we divide Vm and Hm in multiple cores. This will result in more communications, but will make the
#     computation faster when restart is large.
#
#     :param lhs: GlobalMatrix
#     :param rhs: LocallyFullVector
#     :param X0: LocallyFullVector
#     :param restart:
#     :param maxiter:
#     :param tol: relative tolerance.
#     :param atol: absolute tolerance.
#     :param preconditioner: Format: (ID, kwargs (a dict) for the preconditioner)
#     :param COD: Clear Original Data?
#     :param loading_factor: In principle, we will store "loading_factor * 1e8" values in Vm.
#     :return: Return a tuple of 5 outputs:
#
#             1. (LocallyFullVector) results -- The result vector.
#             2. (int) info -- The info which provides convergence information:
#
#                 * 0 : successful exit
#                 * >0 : (convergence to tolerance not achieved) `Info` means number of iterations
#                 * <0 : illegal input or breakdown
#
#             3. (float) beta -- The residual.
#             4. (int) ITER -- The number of outer iterations.
#             5. (str) message
#
#     """
#     assert sIze >= 2, "___mpi_v1_gmres___ routine works for #cores > 2."
#     assert mAster_rank != sEcretary_rank, "This must be the case."
#     assert 0 < loading_factor, f"loading_factor={loading_factor} wrong, it should be in (0, +inf]."
#
#
#     Time_start = MPI.Wtime()
#     A = lhs.M
#     f = rhs.V.toarray().ravel()
#     x0 = np.array(X0.V) # just to make another copy of x0: IMPORTANT TO DO SO!
#     shape0, shape1 = A.shape
#
#     # now we decide how to divide Vm and Hm.
#     Vm_Cor, Vm_Dis, Vm_Ind = ___mpi_v1_distributor___(loading_factor, shape0, restart)
#     # ...
#
#     assert f.shape[0] == x0.shape[0] == shape0 == shape1, "Ax=f shape dis-match."
#
#     if preconditioner is not None:
#         invM = preconditioner.M
#         A = invM @ A
#         f = invM @ f
#
#     if COD: # needs be after preconditioning since the preconditioner will use A.M
#         lhs._M_ = None
#         rhs._V_ = None
#         X0._V_ = None
#
#     ITER = 0
#     BETA = None
#
#     if rAnk == mAster_rank:
#         VV = np.empty((shape0,), dtype=float)
#     else:
#         VV = None
#
#     if rAnk in Vm_Cor:
#         AVJ = np.empty((shape0,), dtype=float)
#         Hm = np.zeros((restart + 1, restart), dtype=float)
#
#         ind = Vm_Cor.index(rAnk)
#         dis = Vm_Dis[ind]
#         Vm = np.empty((dis, shape0), dtype=float)
#
#     else:
#         AVJ = None
#
#     local_range = Vm_Ind[Vm_Cor.index(rAnk)]
#     lrs = local_range.start
#
#     which_cores_need_AVJ = []
#
#     while 1: # always do till break.
#
#         v0 = f - A @ x0
#
#         cOmm.Reduce(v0, VV, op=MPI.SUM, root=mAster_rank)
#
#         if rAnk == mAster_rank:
#             beta = np.sum(VV ** 2) ** 0.5
#             v0 = VV / beta
#         else:
#             beta = None
#
#         cOmm.Bcast([v0, MPI.FLOAT], root=mAster_rank)
#         beta = cOmm.bcast(beta, root=mAster_rank)
#
#         # check stop iteration or not ...
#         if BETA is None: BETA = [beta,] # this is right, do not initial BETA as an empty list.
#         if len(BETA) > 20:
#             BETA = BETA[:1] + BETA[-2:]
#
#         BETA.append(beta)
#         JUDGE, stop_iteration, info = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
#         if stop_iteration: break
#
#         if rAnk == mAster_rank:
#             Vm[0] = v0 # this must be the case, v0 must be stored in the master core.
#         else:
#             v = v0
#
#         for j in range(restart):
#             if j in local_range:
#                 Avj = A @ Vm[j-lrs]
#             else:
#                 Avj = A @ v
#
#
#             cOmm.Reduce(Avj, AVJ, op=MPI.SUM, root=mAster_rank)
#
#
#
#         break
#
#


# def ___mpi_v2_distributor___(restart):
#     """
#
#     :param restart:
#     :return:
#     """
#     # num = restart
#     # parts = sIze
#     # distribution = [num // parts + (1 if x < num % parts else 0) for x in range(parts)][::-1]
#     # print(distribution)
#
#
#     if rAnk== mAster_rank:
#         print(IDict)



def ___mpi_v2_gmres___(lhs, rhs, X0, restart=100, maxiter=20, tol=1e-3, atol=1e-4, preconditioner=None,
                       COD=True, name=None, plot_residuals=False):
    """
    In this version, we divide Vm and Hm in all cores in an optimal way. This will for sure result in more
    communications, but will make the computation faster when restart is large.

    :param lhs: GlobalMatrix
    :param rhs: GlobalVector
    :param X0: LocallyFullVector
    :param restart:
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

    local_ind = list()
    for i in range(restart):
        if (i % sIze) == rAnk:
            local_ind.append(i)

    Time_start = MPI.Wtime()
    A = lhs.M
    f = rhs.V.toarray().ravel()
    x0 = np.array(X0.V) # just to make another copy of x0: IMPORTANT TO DO SO!

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

    if COD: # needs be after preconditioning since the preconditioner will use A.M
        lhs._M_ = None
        rhs._V_ = None
        X0._V_ = None

    ITER = 0
    BETA = None

    AVJ = np.empty((shape0,), dtype=float)
    Hm = np.zeros((restart + 1, restart), dtype=float) # In future, we can change this to sparse matrix.
    if rAnk == mAster_rank:
        VV = np.empty((shape0,), dtype=float)
        Vm = np.empty((restart, shape0), dtype=float)
        HM = np.empty((restart + 1, restart), dtype=float) # use to combine all Hm.
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
        HM = None
        Vm = None
        SUM_Hij_vi = None
        Vs = np.empty((len(local_ind), shape0), dtype=float)

    while 1: # always do till break.

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
                residuals.append(beta)
        JUDGE, stop_iteration, info, JUDGE_explanation = ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...

        if rAnk == mAster_rank:
            Vm[0] = v0
        else:
            Vm = v0
            if 0 in local_ind:
                Vs[local_ind_dict[0]] = v0

        for j in range(restart):

            if rAnk == mAster_rank:
                Avj = A @ Vm[j]
            else:
                Avj = A @ Vm
                if j in local_ind:
                    Vs[local_ind_dict[j]] = Vm

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
                cOmm.Bcast([v_jp1, MPI.FLOAT], root=mAster_rank)
                if rAnk == mAster_rank:
                    Vm[j+1] = v_jp1
                else:
                    Vm = v_jp1

        cOmm.Reduce(Hm, HM, op=MPI.SUM, root=mAster_rank)

        if rAnk == mAster_rank:
            HMT = HM.T
            ls_A = HMT @ HM
            ls_b = HMT[:,0] * beta
            ym = np.linalg.solve(ls_A, ls_b)
            # del HMT, ls_A, ls_b
            x0 += ym.T @ Vm

        cOmm.Bcast([x0, MPI.FLOAT], root=mAster_rank)
        ITER += 1

    x0 = LocallyFullVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    Time_end = MPI.Wtime()

    COST_total = Time_end - Time_start
    message = f" mpi_v2_gmres = [SYSTEM-SHAPE: {A.shape}] [ITER={ITER}][residual=%.2e] costs %.2f, " \
              f"convergence info={info}, restart={restart}, maxiter={maxiter}, " \
              f"stop_judge={JUDGE}: {JUDGE_explanation}]"%(beta, COST_total)

    if rAnk == mAster_rank:
        if plot_residuals:
            ___gmres_plot_residuals___(np.array(residuals), name, 'mpi_gmres_v2')

    return x0, info, beta, ITER, message