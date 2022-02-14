# -*- coding: utf-8 -*-
from SCREWS.frozen import FrozenOnly
# from TOOLS.linear_algebra.solvers.gmres import gmres0 as gmres
import TOOLS.linear_algebra.solvers.serial.scipy_sparse_linalg as spspalinalg_solvers
from TOOLS.linear_algebra.data_structures import GlobalMatrix, DistributedVector, GlobalVector
from time import time
from root.config import *
from scipy import sparse as spspa
from TOOLS.linear_algebra.elementwise_cache import EWC_SparseMatrix
import TOOLS.linear_algebra.__DEPRECATED__.operators as TLO

class _Solver_Base_(FrozenOnly):
    """
    Serve as the parent of all solvers.
    """
    def __init__(self, LS):
        self._LS_ = LS

    def ___PRIVATE_update_outputs___(self, results, messages, info, others):
        """
        Here we update these outputs:

            1. results -- The result vector
            2. (str) messages -- A message (which normally summarize the solving process).
            3. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * <0 : illegal input or breakdown

            4. others -- All other things

        :param results:
        :param messages:
        :param int info:
        :param others:
        :return:
        """
        assert results.__class__.__name__ == 'DistributedVector'
        self._LS_.solve._results_, \
        self._LS_.solve._message_, \
        self._LS_.solve._info_, \
        self._LS_.solve._others_ = results, messages, info, others

    def ___PRIVATE_check_regularity_compatibility___(self, allowed):
        """
        Check if the regularity of current linear system is allowed.

        :param allowed: The pool of all allowed regularities.
        :return:
        """
        arr = list()
        for ar in allowed:
            ar = ar.replace(' ', '')
            ar = ar.replace('\n', '')
            arr.append(ar)
        if self._LS_.structure.regularity_ravel in arr:
            return arr.index(self._LS_.structure.regularity_ravel)
        else:
            if rAnk == mAster_rank:
                error_str = f"Linear system <{self._LS_.name}>'s regularity:\n\n" \
                            f"{self._LS_.structure.regularity}\n" \
                            f"does not match one of allowed ones:\n"
                arr = list()
                for ar in allowed:
                    ar = ar.replace('        ', '')
                    arr.append(ar)
                for ar in arr:
                    error_str += ar
                raise Exception(error_str)
            else:
                raise Exception()





class Direct(_Solver_Base_):
    """``GMRES`` solver is a general solver using the GMRES scheme."""
    def __init__(self, LS):
        super(Direct, self).__init__(LS)
        self._freeze_self_()

    @property
    def _default_routine_(self):
        return 'spspalinalg'

    def spspalinalg(self):
        """"""
        ts = time()
        assert self._LS_.structure.IS_all_global_matrices, \
            "GMRES-MPI1 routine only works for IS_all_global_matrices linear system."
        A = self._LS_.lhs.___PRIVATE_convert_into_one_sparse_A___()
        AMS = A.M.shape
        b = self._LS_.rhs.___PRIVATE_convert_into_one_sparse_b___()

        A = A.___PRIVATE_gather_M_to_core___(clean_local=True)
        A = GlobalMatrix(A)

        b = b.___PRIVATE_gather_V_to_core___(clean_local=True)
        b = GlobalVector(b)

        results, info, TOL, ITER = spspalinalg_solvers.spsolve(A, b)[:4]
        messages = f"[gmres] A shape={AMS}, Solver: scipy.sparse.linalg.spsolve, " \
                   f"time cost={int((time()-ts)*100)/100} seconds."
        others = None
        self.___PRIVATE_update_outputs___(results, messages, info, others)





class ITERATIVE(_Solver_Base_):
    def __init__(self, LS):
        super(ITERATIVE, self).__init__(LS)
        self._freeze_self_()

    @property
    def _default_routine_(self):
        return 'MPI1'

    def spspalinalg(self, x0, restart=100, maxiter=1000, tol=1e-5, solver='gcrotmk', conserving=False):
        """
        We use the gmres scheme to do the solving.

        It is the first gmres solver and it not optimal at efficiency. So it is called MPI1.

        :param x0:
        :param restart:
        :param maxiter:
        :param tol:
        :param solver:
        :param conserving:

        """
        ts = time()
        assert self._LS_.structure.IS_all_global_matrices, \
            "GMRES-MPI1 routine only works for IS_all_global_matrices linear system."
        A = self._LS_.lhs.___PRIVATE_convert_into_one_sparse_A___()
        AMS = A.M.shape
        b = self._LS_.rhs.___PRIVATE_convert_into_one_sparse_b___()

        A.___PRIVATE_gather_M_to_core___(clean_local=True)
        b._V_ = b.___PRIVATE_gather_V_to_core___(clean_local=True)


        if conserving or solver == 'spsolve':

            results, info, TOL, ITER = spspalinalg_solvers.spsolve(A, b)[:4]

        else:

            assert x0.__class__.__name__ in ('DistributedVector',)

            if solver  == 'gmres':
                results, info, TOL, ITER = spspalinalg_solvers.gmres(
                    A, b, x0, restart=restart, maxiter=maxiter, tol=tol)[:4]
            elif solver == 'gcrotmk':
                results, info, TOL, ITER = spspalinalg_solvers.gcrotmk(
                    A, b, x0, maxiter=maxiter, tol=tol)[:4]
            elif solver == 'bicgstab':
                results, info, TOL, ITER = spspalinalg_solvers.bicgstab(
                    A, b, x0, maxiter=maxiter, tol=tol)[:4]
            else:
                raise NotImplementedError()

        messages = f"[{solver}] A shape={AMS}, Solver: GMRES, routine: MPI1, " \
                   f"convergence info={info}, residual={TOL}, restart={restart}, " \
                   f"maxiter={maxiter}, used iterations={ITER}, " \
                   f"time cost={int((time()-ts)*100)/100} seconds."
        others = None
        self.___PRIVATE_update_outputs___(results, messages, info, others)




class icpsNS(_Solver_Base_):
    """
    ``GMRES`` solver is a general solver using the GMRES scheme.
    """
    def __init__(self, LS):
        super(icpsNS, self).__init__(LS)
        self.DO_reset_cache()
        self.___PRIVATE_OS1_reset_cache___()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._OS1_GMA_ = None

    @property
    def _default_routine_(self):
        raise Exception("Has to provide a routine name.")

    def ___PRIVATE_OS1_reset_cache___(self):
        self._OS1_A00_ = None
        self._OS1_A01_ = None
        self._OS1_A10_ = None
        self._OS1_A11_ = None

    def ___PRIVATE_OS1_A_generator___(self, i):
        """"""
        A00 = self._OS1_A00_[i]
        A01 = self._OS1_A01_[i]
        A10 = self._OS1_A10_[i]
        A11 = self._OS1_A11_
        return spspa.bmat(([A00, A01],
                           [A10, A11]), format='csc')

    def OS1(self, x0, restart=100, maxiter=1000, tol=1e-5, solver='gcrotmk', conserving=False):
        """

        :param x0:
        :param restart:
        :param maxiter:
        :param tol:
        :param solver:
        :param conserving:
        :return:
        """
        ts = time()
        allowed_regularities = (
        """
        EWC_SpaseMat EWC_SpaseMat -----GM----- -----GM----- | EWC_ColVec
        EWC_SpaseMat ------------ ------------ ------------ | gv
        -----GM----- ------------ -----GM----- ------------ | gv
        -----GM----- ------------ ------------ -----GM----- | gv
        """
        ,)
        FORMAT = self.___PRIVATE_check_regularity_compatibility___(allowed_regularities)
        if FORMAT not in (0,): raise NotImplementedError()

        #  |alpha  beta|  |g|
        #  |gamma  F   |  |h|
        #  alpha = lhs(0:2, 0:2)
        A00 = self._LS_.lhs[0][0]
        A01 = self._LS_.lhs[0][1]
        A10 = self._LS_.lhs[1][0]
        assert A00.elements == A01.elements == A10.elements == self._LS_.rhs[0].elements
        if rAnk == mAster_rank:
            l1, r1 = np.shape(A01[0])[1], np.shape(A10[0])[0]
        else:
            l1, r1 = None, None
        l1, r1 = cOmm.bcast([l1, r1], root=mAster_rank)
        Z11 = spspa.csc_matrix((r1, l1))
        self._OS1_A00_ = A00
        self._OS1_A01_ = A01
        self._OS1_A10_ = A10
        self._OS1_A11_ = Z11

        alpha_inv = EWC_SparseMatrix(A00.elements, self.___PRIVATE_OS1_A_generator___, 'no_cache')
        alpha_inv = alpha_inv.inv

        if self._OS1_GMA_ is None:
            GMP3, GMu2 = A10.gathering_matrices
            self._OS1_GMA_ = GMu2.GMs[0].DO_hstack(GMP3.GMs[0])

        alpha_inv.gathering_matrices = (self._OS1_GMA_, self._OS1_GMA_)
        alpha_inv = alpha_inv.assembled

        len_1 = A10.gathering_matrices[0].GLOBAL_num_dofs
        len_2 = self._LS_.lhs[0][2].shape[1]
        beta  = TLO.bmat(([self._LS_.lhs[0][2]         , self._LS_.lhs[0][3]],
                          [GlobalMatrix((len_1, len_2)), None               ]), format='csc')
        gamma = TLO.bmat(([self._LS_.lhs[2][0], GlobalMatrix((len_2, len_1))],
                          [self._LS_.lhs[3][0], None                        ]), format='csc')
        F     = TLO.bmat(([self._LS_.lhs[2][2], None                        ],
                          [None               , self._LS_.lhs[3][3]         ]), format='csc')
        B0 = self._LS_.rhs[0].assembled
        g = TLO.concatenate((B0              ,self._LS_.rhs[1]))
        h = TLO.concatenate((self._LS_.rhs[2],self._LS_.rhs[3]))

        if saFe_mode:
            assert alpha_inv.shape[0] == beta.shape[0]  == g.shape[0], "[icpsNS] -> [OS1]"
            assert gamma.shape[0]     == F.shape[0]     == h.shape[0], "[icpsNS] -> [OS1]"
            assert alpha_inv.shape[1] == gamma.shape[1] == g.shape[0], "[icpsNS] -> [OS1]"
            assert beta.shape[1]      == F.shape[1]     == h.shape[0], "[icpsNS] -> [OS1]"

        alpha_inv._M_ = alpha_inv.M.tocsr()
        alpha_inv.IS_regularly_distributed = 'row'

        # make the operations done in single core ...
        alpha_inv.___PRIVATE_gather_M_to_core___(clean_local=True)
        beta.___PRIVATE_gather_M_to_core___(clean_local=True)
        gamma.___PRIVATE_gather_M_to_core___(clean_local=True)
        F.___PRIVATE_gather_M_to_core___(clean_local=True)
        # ...

        A1 = alpha_inv @ beta
        A2 = gamma @ A1
        del A1
        lhs = F - A2
        del A2
        rhs = h - gamma @ (alpha_inv @ g)

        cOmm.barrier()
        total_shape = g.shape[0] + h.shape[0]
        del F, gamma, h
        non_zeros = lhs.GLOBAL_approx_nnz
        lhs_shape = lhs.shape
        total_entries = lhs_shape[0] * lhs_shape[1]
        assert lhs_shape[0] <= non_zeros <= total_entries, "Something wrong."
        sparsity = (1 - non_zeros / total_entries) * 100
        messages = f'[icpsNS-OS1] Total system shape {(total_shape, total_shape)}, ' \
                   f'Reduced system shape ' + str(lhs_shape) + \
                    ' @sparsity~={}% with ~ {} non-zeros.\n'.format('%0.2f'%sparsity, non_zeros)

        t1 = time()
        if conserving or solver  == 'spsolve':
            res_q, info, TOL, ITER = spspalinalg_solvers.spsolve(lhs, rhs)[:4]
        else:
            if solver  == 'gmres':
                res_q, info, TOL, ITER = spspalinalg_solvers.gmres(
                    lhs, rhs, x0, restart=restart, maxiter=maxiter, tol=tol)[:4]
            elif solver == 'gcrotmk':
                res_q, info, TOL, ITER = spspalinalg_solvers.gcrotmk(
                    lhs, rhs, x0, maxiter=maxiter, tol=tol)[:4]
            elif solver == 'bicgstab':
                res_q, info, TOL, ITER = spspalinalg_solvers.bicgstab(
                    lhs, rhs, x0, maxiter=maxiter, tol=tol)[:4]
            else:
                raise NotImplementedError()

        reduced_system_solving_time = time()- t1
        del lhs, rhs

        res_up = alpha_inv @ (g - beta @ res_q)
        del alpha_inv, beta

        messages += 'Solver costs {}s among which {}s is for solving the reduced system.\n'.format(
                    '%0.2f'%(time() - ts), '%0.2f'%reduced_system_solving_time) + \
                   f"[{solver}] on the reduced system: " + \
                   f"convergence info={info}, residual={TOL}, " \
                   f"restart={restart}, maxiter={maxiter}, used iterations={ITER}"
        others = [res_q,] # Lets save this as x0 for next iteration.
        res_up = res_up.___PRIVATE_gather_V_to_core___()
        if rAnk == mAster_rank:
            res_up = spspa.csc_matrix(res_up[:,np.newaxis])
        else:
            res_up = spspa.csc_matrix((g.shape[0], 1))
        res_up = DistributedVector(res_up)
        results = TLO.concatenate((res_up, res_q))
        self.___PRIVATE_update_outputs___(results, messages, info, others)