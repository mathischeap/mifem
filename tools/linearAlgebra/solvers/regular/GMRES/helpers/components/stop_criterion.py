# -*- coding: utf-8 -*-


def ___gmres_stop_criterion___(tol, atol, ITER, maxiter, BETA):
    """
    :param tol: relative tolerance.
    :param atol: absolute tolerance
    :param ITER:
    :param maxiter: int, str
        A positive integer.

        if maxiter is a str, it must be a numeric str, and it means it is a
        strong maxiter, that is no matter what happened, we will iterate the
        solver for this many times. So it is a forced amount of iterations.

    :param BETA: A list of beta (residual) of some recent iterations.
    :return:
    """

    if isinstance(maxiter, str):

        judge_1 = False
        judge_3 = False
        judge_4 = False
        judge_5 = False

        MAXITER = int(maxiter)
        judge_2 = ITER >= MAXITER # judge 2: reach max iteration number


    else:

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
        if len(BETA) >= 5:

            if len(BETA) > 10:
                CHECKS = BETA[-8:]
            else:
                CHECKS = BETA[1:]

            JUDGES = list()

            Lc = len(CHECKS)

            for i in range(Lc-1):
                j = i + 1

                OLD = CHECKS[i]
                NEW = CHECKS[j]

                progress = OLD - NEW

                if progress < 0 :
                    JUDGES.append(True)

                else:
                    if progress / OLD < tol: # slow converging
                        JUDGES.append(True)
                    else:
                        JUDGES.append(False)

            cT = JUDGES.count(True)
            cF = JUDGES.count(False)

            if cT >= 3:
                judge_5 = True
            elif cT >= cF:
                judge_5 = True
            else:
                judge_5 = False

        else:
            judge_5 = False

    # ................................................................

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