# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/14/2022 9:28 PM
"""


def nLS_stop_criterion(BETA, atol, ITER, maxiter):
    """

    Parameters
    ----------
    BETA
    atol
    ITER
    maxiter : int, str
        If it is a str, then it is a forced maxiter, which is the only criterion of stopping
        the iterations.

    Returns
    -------

    """

    if isinstance(maxiter, str):

        MAXITER = int(maxiter)
        judge_2 = ITER >= MAXITER  # judge 2: reach max iteration number

        judge_1 = False
        judge_3 = False
        judge_4 = False

    else:

        beta = BETA[-1]
        judge_1 = beta < atol  # judge 1: reach absolute tolerance.
        judge_2 = ITER >= maxiter  # judge 2: reach max iteration number

        # judge 3: divergence
        if len(BETA) > 1 and BETA[-1] > BETA[-2]:  # error grows after one iteration
            if BETA[-2] > 1 and (BETA[-1] - BETA[-2]) > 0.5 * BETA[-2]:
                judge_3 = True
            elif BETA[-1] > 10e6:
                judge_3 = True
            elif (BETA[-1] - BETA[-2]) > 10:
                judge_3 = True
            elif BETA[-1] > 10 * BETA[-2]:
                judge_3 = True
            else:
                judge_3 = False
        else:
            judge_3 = False

        # judge_4: slow converging
        judge_4 = False  # we currently turn off this judge, if it is slow, we just do more iterations.

        # if len(BETA) > 1 and beta < BETA[-2]:
        #     beta_old = BETA[-2]
        #     progress = beta_old - beta
        #     if progress / beta_old < 0.001:  # slow converging
        #         judge_4 = True
        #     else:
        #         judge_4 = False
        # else:
        #     judge_4 = False

    # -------------------------------------------------------------------------------
    if judge_1 or judge_2 or judge_3 or judge_4:

        stop_iteration = True

        if judge_1:  # reach atol
            info = 0
            JUDGE = 1
            JUDGE_explanation = 'reach absolute tol'

        elif judge_2:  # reach maxiter
            info = ITER
            JUDGE = 2
            JUDGE_explanation = 'reach maxiter'

        elif judge_3:  # diverging
            info = -1
            JUDGE = 3
            JUDGE_explanation = 'diverging'

        elif judge_4:  # very slow converging; the progress is lower than the tol
            info = ITER
            JUDGE = 4
            JUDGE_explanation = 'very slow converging'
        else:
            raise Exception()

    else:  # do not stop iterations.
        stop_iteration = False
        info = None
        JUDGE = 0
        JUDGE_explanation = ''

    assert stop_iteration in (True, False), "stop_iteration has to be set."

    return JUDGE, stop_iteration, info, JUDGE_explanation