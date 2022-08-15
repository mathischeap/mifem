# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/4 21:19
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.solvers.regular.allocator import RegularSolverDistributor
from tools.linear_algebra.nonlinear_system.solve.Newton_Raphson.helpers.stop_criterion import nLS_stop_criterion

from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate
from time import time

class nLS_Solve_NR_regular(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._freeze_self_()

    def __call__(self, x0, atol=1e-4, maxiter=10,
        LS_solver_para='GMRES', LS_solver_kwargs=None):
        """

        Parameters
        ----------
        x0 :
            The initial guess
        atol :
            Absolutely tolerance; norm of delta_x.
        maxiter : int
            A positive integer.
        LS_solver_para : str, tuple
            The parameters to select and config the linear solver.

            We will use `LS_solver_para` to find the correct allocator and then parse the
            right solver.

            A tuple like:
                ('GMRES', {routine='auto', name="A-FANCY-solver"})

                We will find the way to the regular linear solver allocator and use
                'GMRES', routine='auto', name="A-FANCY-solver"
                as the (kw)args. This will give a linear GMRES solver named "A-FANCY-solver"
                which will automatically choose routine.

        LS_solver_kwargs : dict, NoneType
            The kwargs sent to the linear solver.
            (If needed, x0 will be sent as for the first Newton-Raphson iteration.)

       Returns
       -------

       """
        t_start = time()
        assert isinstance(atol, (int, float)) and atol > 0, f"atol={atol} is wrong."
        if LS_solver_kwargs is None: LS_solver_kwargs = dict()
        LS_solver_kwargs['COD'] = False # make sure we do not automatically clear data.

        if isinstance(LS_solver_para, str): LS_solver_para = (LS_solver_para, dict())
        assert isinstance(LS_solver_para, tuple), f"please put LS_solver_para in a tuple."
        assert isinstance(LS_solver_kwargs, dict), f"please put LS_solver_kwargs in a dict."

        linear_solver_name = LS_solver_para[0]
        linear_solver_parameters = LS_solver_para[1]

        #--find the linear solver caller; call it with LS_solver_kwargs to solve a linear system--
        if linear_solver_name in RegularSolverDistributor.___solver_name___():
            linear_solver_caller = RegularSolverDistributor(linear_solver_name,
                                                            **linear_solver_parameters)
        else:
            raise NotImplementedError()

        #----- initialize the variables -------------------------------------------------------1

        xi = x0
        ITER = 0
        message = ''
        BETA = list()

        #------------ define intermediate variables -------------------------------------------1

        itmV = list()
        for uv in self._nLS_.unknown_variables:
            itmV.append(uv.shadow)

        t_iteration_start = time()
        while 1:
            ITER += 1

            #---------- distribute xi to the intermediate forms -------------------------------2
            xi.do.distributed_to(*itmV)

            #--- evaluate blocks of the left-hand-side matrix with xi -------------------------2
            S0, S1 = self._nLS_.shape
            LHS = [[list() for _ in range(S1)] for _ in range(S0)]
            for i in range(S0):

                NTi = self._nLS_.nonlinear_terms[i]
                tv = self._nLS_.test_variables[i]

                for j in range(S1):
                    #__________ we first find the linear term from A _________________________3
                    linear_term = self._nLS_.A[i][j]
                    if linear_term is not None:
                        LHS[i][j].append(linear_term)
                    else:
                        pass

                    #____ we then look at the nonlinear terms to find the contributions _____3

                    uv = self._nLS_.unknown_variables[j]

                    if NTi == 0 or NTi is None:
                        pass

                    else:

                        known_pairs = list()
                        for _ in self._nLS_.unknown_variables:
                            if _ is not uv:
                                known_pairs.append((_, self.___PRIVATE_find_shadow___(_, itmV)))

                        uv_pair = (uv, self.___PRIVATE_find_shadow___(uv, itmV))

                        if isinstance(NTi, (list, tuple)):
                            for nti in NTi:
                                contribution_2_Aij = \
                                    nti.do.derivative_contribution_to(tv, uv_pair, *known_pairs)
                                if contribution_2_Aij is not None:
                                    LHS[i][j].append(contribution_2_Aij)
                        else:
                            contribution_2_Aij = \
                                NTi.do.derivative_contribution_to(tv, uv_pair, *known_pairs)
                            if contribution_2_Aij is not None:
                                LHS[i][j].append(contribution_2_Aij)

                    #----------- sum up blocks ____________________________________________3
                    LEN = len(LHS[i][j])
                    if LEN == 0:
                        # noinspection PyTypeChecker
                        LHS[i][j] = None
                    elif LEN == 1:
                        # noinspection PyTypeChecker
                        LHS[i][j] = LHS[i][j][0]
                    elif LEN == 2:
                        # noinspection PyTypeChecker
                        LHS[i][j] = LHS[i][j][0] + LHS[i][j][1]
                    else:
                        # noinspection PyTypeChecker
                        LHSij0 = LHS[i][j][0]
                        LHS[i][j] = LHSij0.___PRIVATE_sum___(LHS[i][j][1:])

            #--------- now we evaluate vector f ---------------------------------------------2
            f = self._nLS_.do.evaluate_f(itmV)

            #--------- A, f -----------------------------------------------------------------2
            A = bmat(LHS)
            f = - concatenate(f) # notice the minus here, it is important.

            A.gathering_matrices = (self._nLS_.test_variables, self._nLS_.unknown_variables)
            f.gathering_matrix = self._nLS_.test_variables

            #---------------- adopt customizations ------------------------------------------2
            customizations = self._nLS_.customize.customizations

            for customization in customizations:
                method_name = customization[0]
                method_para = customization[1]

                if method_name == "set_no_evaluation":

                    A.customize.identify_global_row(method_para)
                    f.customize.set_assembled_V_i_to(-1, 0)

                else:
                    raise NotImplementedError(f"Cannot handle customization = {method_name}")

            #----------- assemble ---------------------------------------------------------2

            A = A.assembled
            f = f.assembled

            #---------- solve -------------------------------------------------------------2
            if linear_solver_caller._solver_.__class__.__name__ == 'Direct':

                LSR = linear_solver_caller(A, f, **LS_solver_kwargs)

            else:

                LSR = linear_solver_caller(A, f, 0, **LS_solver_kwargs)

            dx, LSm = LSR[0], LSR[4]

            #------------ check the progress ----------------------------------------------2
            BETA.append(dx.vector_norm)
            if len(BETA) > 10: BETA = [BETA[0],] + BETA[-2:]
            JUDGE, stop_iteration, convergence_info, JUDGE_explanation = \
                nLS_stop_criterion(BETA, atol, ITER, maxiter)

            # ---------- here we use known solution xi to compute dx (xi1 = xi + dx) ------2
            xi1 = xi + dx

            if stop_iteration:
                break
            else:
                xi = xi1

        #----------------------------------------------------------------------------------------1
        t_iteration_end = time()
        Ta = t_iteration_end-t_start
        TiT = t_iteration_end-t_iteration_start

        #----------------------------------------------------------------------------------------1
        message += f"<nonlinear_solver>" \
                   f" = [RegularNewtonRaphson-{A.shape}]" \
                   f" of {self._nLS_.nonlinear_terms.num} nonlinear terms" \
                   f" = Linear solver args: {LS_solver_para} {LS_solver_kwargs}" \
                   f" = [ITER:{ITER}]" \
                   f" = [beta-{BETA[-1]}]" \
                   f" = [{convergence_info}-{JUDGE_explanation}]" \
                   f" >>> nLS solving costs %.2f, each ITER cost %.2f"%(Ta, TiT/ITER)\
                   + '\n >>> Last Linear Solver Message:\n' + LSm

        #=======================================================================================1

        return xi1, convergence_info, BETA[-1], ITER, message


    def ___PRIVATE_find_shadow___(self, obj, shadows):
        """"""
        index = self._nLS_.unknown_variables.index(obj)
        return shadows[index]





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
