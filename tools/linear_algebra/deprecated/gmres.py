# -*- coding: utf-8 -*-
from root.config import *
from scipy import sparse as spspa
from scipy.sparse import linalg as spspalinalg
from tools.linear_algebra.data_structures.global_matrix.main import DistributedVector
from screws.frozen import FrozenOnly
from screws.exceptions import LinerSystemSolverDivergenceError



def gmres0(AA, bb, X0, restart=100, maxiter=1000, tol=1e-4):
    """
    The first gmres scheme.

    :param GlobalMatrix AA:
    :param GlobalVector bb:
    :param DistributedVector X0:
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
    """
    assert AA.__class__.__name__ == 'GlobalMatrix'
    assert bb.__class__.__name__ == 'GlobalVector'
    assert X0.__class__.__name__ == 'DistributedVector'
    assert maxiter >= 1, "maxiter must be >= 1."
    assert restart >= 3, "restart must be >= 3."
    assert tol > 0, "tol must be > 0."

    LNC = AA.nonempty_columns # the columns have non-zero values.
    ANC = cOmm.gather(LNC, root=mAster_rank) # all stored in master core.
    distribute_vector = ___gmres_Distribute_Vector___(AA, ANC, LNC)
    combine_vector = ___gmres_Combine_Vector___()
    del ANC, LNC

    A = AA.M
    f = bb.V
    x0 = X0.V

    ITER = 0
    BETA = None
    while 1:
        r0 = f - A @ x0
        beta = ___gmres_norm___(r0)

        # if rAnk == 0:
        #     print('gmres0:', ITER, 'restart:', restart, 'error:', beta, flush=True)


        # check stop iteration or not ...
        if BETA is None: BETA = [beta,]
        if len(BETA) > 20: BETA = BETA[-5:]
        BETA.append(beta)
        stop_iteration, info = ___gmres_stop_criterion___(tol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...


        v0 = r0 / beta
        v0 = combine_vector(v0) # v0 only in master core now.
        if rAnk == mAster_rank:
            Vm = np.empty((restart, v0.shape[0]))
            Vm[0] = v0
            Hm = spspa.lil_matrix((restart+1, restart))
        else:
            pass

        if rAnk != mAster_rank: vj = None

        for j in range(restart):
            if rAnk == mAster_rank:
                vj = Vm[j]
            vj = distribute_vector(vj)
            wj = A @ vj
            wj = combine_vector(wj)
            if rAnk == mAster_rank:
                hij_vi = None
                for i in range(0, j+1):
                    Hm_ij = np.sum(wj * Vm[i])
                    if hij_vi is None:
                        # noinspection PyUnresolvedReferences
                        hij_vi = Hm_ij * Vm[i]
                    else:
                        hij_vi += Hm_ij * Vm[i]
                    Hm[i, j] = Hm_ij
                hat_v_jp1 = wj - hij_vi
                Hm[j+1, j] = np.sum(hat_v_jp1**2) ** 0.5
                if j < restart-1:
                    Vm[j+1] = hat_v_jp1 / Hm[j+1, j]
            else:
                pass

        if rAnk == mAster_rank:
            Hm = Hm.tocsr()
            HmT = Hm.T
            ls_A = HmT @ Hm
            ls_b = HmT[:,0] * beta
            ym = spspalinalg.spsolve(ls_A, ls_b)

            del HmT, ls_A, ls_b
            Vm_ym = Vm.T @ ym
        else:
            Vm_ym = None

        Vm_ym = distribute_vector(Vm_ym)
        x0 += Vm_ym
        ITER += 1

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres0 diverges after {ITER} iterations with error reaching {beta}.")

    return x0, info, beta, ITER

def ___gmres_stop_criterion___(tol, ITER, maxiter, BETA):
    """
    :param tol:
    :param ITER:
    :param maxiter:
    :param BETA:
    :return:
    """
    beta = BETA[-1]
    # judge 1: reach tol
    judge_1 = beta < tol
    # judge 2: reach max iteration number
    judge_2 = ITER >= maxiter
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
    # ...
    if judge_1 or judge_2 or judge_3:
        if judge_1: # reach tol
            info = 0
        elif judge_2: # reach maxiter
            info = ITER
        elif judge_3: # divergence
            info = -1
        stop_iteration = True
    else:
        info = None
        stop_iteration = False
    return stop_iteration, info

def ___gmres_norm___(r):
    r = cOmm.gather(r, root=mAster_rank)
    if rAnk == mAster_rank:
        r = np.sum(r)
        # noinspection PyUnresolvedReferences
        norm = np.sum(r.power(2))**0.5
    else:
        norm = None
    norm = cOmm.bcast(norm, root=mAster_rank)
    return norm

class ___gmres_Combine_Vector___(FrozenOnly):
    def __init__(self):
        self._tree_ = [Hi for Hi in tRee(2)]
        self._freeze_self_()

    def __call__(self, v):
        for Hi in self._tree_:
            if Hi is None:
                pass
            elif Hi[0] == 'send':
                cOmm.send(v, **Hi[1])
                v = None
            elif Hi[0] == 'recv':
                v += cOmm.recv(**Hi[1])
            else:
                raise Exception()
        if rAnk == mAster_rank: v = v.toarray()[:,0]
        return v

class ___gmres_Distribute_Vector___(FrozenOnly):
    def __init__(self, AA, ANC, LNC):
        self._ANC_ = ANC
        self._shape_ = (AA.shape[0], 1)
        self._LI_ = LNC
        self._indptr_ = [0, len(self._LI_)]
        self._freeze_self_()

    def __call__(self, v):
        if rAnk == mAster_rank:
            V = [v[nnc] for nnc in self._ANC_]
        else:
            V = None
        V = cOmm.scatter(V, root=mAster_rank)
        V = spspa.csc_matrix((V, self._LI_, self._indptr_), shape=self._shape_)
        return V






def ___DEPRECATED_gmres_combine_vector___(v):
    """
    We gather v from all cores into master and sum them up.

    :param v:
    :return:
    """
    v = cOmm.gather(v, root=mAster_rank)
    if rAnk == mAster_rank:
        # noinspection PyUnresolvedReferences
        v = np.sum(v).toarray()[:,0]
    return v

def ___DEPRECATED_gmres_distribute_vector___(v, nonempty_columns, ndofs, indices):
    """
    Distribute vector to all cores.

    :param v:
    :return:
    """
    # vi = None
    # if rAnk == mAster_rank:
    #     for i, nnc in enumerate(nonempty_columns):
    #         if v.__class__.__name__ == 'csc_matrix':
    #             VI = spspa.csc_matrix((v[nnc, 0].T.toarray()[0], indices[nnc], [0, np.count_nonzero(nnc)]),
    #                               shape=(ndofs, 1))
    #         else:
    #             VI = spspa.csc_matrix((v[nnc], indices[nnc], [0, np.count_nonzero(nnc)]),
    #                               shape=(ndofs, 1))
    #         if i != mAster_rank:
    #             cOmm.send(VI, dest=i, tag=i)
    #         else:
    #             vi = VI
    # else:
    #     vi = cOmm.recv(source=mAster_rank, tag=rAnk)
    # return vi
    if rAnk == mAster_rank:
        V = list()
        for i, nnc in enumerate(nonempty_columns):
            VI = spspa.csc_matrix((v[nnc], indices[nnc], [0, np.count_nonzero(nnc)]), shape=(ndofs, 1))
            V.append(VI)
    else:
        V = None
    V = cOmm.scatter(V, root=mAster_rank)
    return V







def gmres1(AA, bb, X0, restart=100, maxiter=1000, tol=1e-4):
    """
    The second gmres scheme (hope it is better than the first one).

    This idea for this solver is that we break vector vj and store it in each cores. Therefore,
    we can avoid collecting and distributing vector vj. But, to do the A @ vj, we have to bcast
    vj in each core to all cores which may be slow. But overall, I believe this is faster.

    It turns out to be very very slow! do NOT use this routine! I over-estimated the speed of
    communication between cores! It is not that fast.

    :param GlobalMatrix AA:
    :param GlobalVector bb:
    :param DistributedVector X0:
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
    """
    assert AA.__class__.__name__ == 'GlobalMatrix'
    assert bb.__class__.__name__ == 'GlobalVector'
    assert X0.__class__.__name__ == 'DistributedVector'

    assert maxiter >= 1, "maxiter must be >= 1."
    assert restart >= 3, "restart must be >= 3."
    assert tol > 0, "tol must be > 0."

    # if AA.mtype != 'csc': AA.do.tocsc()

    A = AA.M
    f = bb.V
    x0 = X0.V

    # nonempty_columns = np.diff(A.indptr) != 0 # the columns have non-zero values.
    # nonempty_columns = cOmm.gather(nonempty_columns, root=mAster_rank) # all stored in master core.
    # if rAnk == mAster_rank:
    #     _, ndofs = np.shape(nonempty_columns)
    #     assert _ == sIze
    #     indices = np.arange(0, ndofs)
    # else:
    #     ndofs = indices = None

    LEN = A.shape[0]
    assert LEN == A.shape[1] == f.shape[0] == x0.shape[0], "Shape dis-match."
    # vector distribution
    allocated = [LEN // sIze + (1 if x < LEN % sIze else 0) for x in range(sIze)]
    LVR = np.arange(sum(allocated[0:rAnk]), sum(allocated[:(rAnk+1)]))
    group_num = np.sqrt(sIze)
    if group_num % 1 > 0.75:
        group_num = int(group_num) + 1
    else:
        group_num = int(group_num)
    CG = gRoup_cores(-1, group_num=group_num)
    AVR = list()
    for i in range(sIze): AVR.append(range(sum(allocated[:i]), sum(allocated[:(i + 1)])))

    # print(rAnk, CG)

    if rAnk in CG['leaders']: # is a group leader
        GVR = [range(sum(allocated[:rAnk]), sum(allocated[:(rAnk + 1)])),]
        for i in CG[rAnk]:
            GVR.append(range(sum(allocated[:i]), sum(allocated[:(i + 1)])))
        # print(GVR)
    else:
        GVR = None

    # will be abandoned ...
    # if rAnk == mAster_rank:
    #     AllocateD = allocated
    #     Vec_rangE = list()
    #     for i in range(sIze):
    #         Vec_rangE.append(np.arange(sum(allocated[:i]), sum(allocated[:(i + 1)])))
    # else:
    #     AllocateD, Vec_rangE = None, None
    # ...

    GD = _gmres1_GROUP_gather_and_redistribute_vector(CG, LVR, AVR, GVR, LEN)

    ITER = 0
    BETA = None

    f = GD(f)
    # f = _gmres1_gather_and_redistribute_vector(f, Vec_rangE, AllocateD)
    x0 = GD(x0)
    # x0 = _gmres1_gather_and_redistribute_vector(x0, Vec_rangE, AllocateD)
    while 1:
        r0 = f - _gmres1_MAT_dot_VEC(A, x0)
        r0 = GD(r0)
        beta = _gmres1_norm_of_vector(r0)

        # if rAnk == 0:
        #     print('gmres1:', ITER, 'restart:', restart, 'error:', beta, flush=True)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,]
        if len(BETA) > 20: BETA = BETA[-5:]
        BETA.append(beta)
        stop_iteration, info = ___gmres_stop_criterion___(tol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...
        v0 = r0 / beta
        # v0 = GD(v0)
        # v0 = _gmres1_gather_and_redistribute_vector(v0, Vec_rangE, AllocateD)
        if rAnk == mAster_rank:
            Hm = spspa.lil_matrix((restart+1, restart))
        Vm = [v0,]


        for j in range(restart):
            vj = Vm[j]
            wj = _gmres1_MAT_dot_VEC(A, vj)
            wj = GD(wj)
            # wj = _gmres1_gather_and_redistribute_vector(wj, Vec_rangE, AllocateD)
            hij_vi = None
            for i in range(0, j + 1):
                Hm_ij = np.sum(wj.multiply(Vm[i]))
                Hm_ij = cOmm.allreduce(Hm_ij, op=MPI.SUM)
                if rAnk == mAster_rank: Hm[i,j] = Hm_ij
                if hij_vi is None:
                    # noinspection PyUnresolvedReferences
                    hij_vi = Hm_ij * Vm[i]
                else:
                    hij_vi += Hm_ij * Vm[i]
            hat_v_jp1 = wj - hij_vi
            Hm_j1_j = _gmres1_norm_of_vector(hat_v_jp1)
            if rAnk == mAster_rank: Hm[j+1,j] = Hm_j1_j
            if j < restart-1:
                Vm.append(hat_v_jp1 / Hm_j1_j)

        if rAnk == mAster_rank:
            # print(Hm.toarray())
            Hm = Hm.tocsr()
            HmT = Hm.T
            ls_A = HmT @ Hm
            ls_b = HmT[:,0] * beta
            ym = spspalinalg.spsolve(ls_A, ls_b)
            del HmT, ls_A, ls_b
        else:
            ym = np.empty(restart, dtype='d')
        cOmm.Bcast([ym, MPI.DOUBLE], root=mAster_rank)
        Vm = spspa.hstack(Vm)
        x0 += spspa.csc_matrix((Vm @ ym)[:, np.newaxis])
        ITER += 1

    x0 = DistributedVector(x0)

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres1 diverges after {ITER} iterations with error reaching {beta}.")

    return x0, info, beta, ITER


def _gmres1_MAT_dot_VEC(M, v):
    """
    w = M @ v.

    :param scipy.sparse.csc_matrix M: Randomly distributed sparse matrix.
    :param v: Redistributed vector.
    :return:
    """
    w = None
    for i in range(sIze):
        lv = cOmm.bcast(v, root=i)
        if w is None:
            w = M @ lv
        else:
            w += M @ lv
    # v = cOmm.allgather(v)
    # v = np.sum(v)
    # w = M @ v
    return w

def _gmres1_norm_of_vector(v):
    """
    Compute the norm of a vector.

    :param v: The vector to compute the norm.
    :return:
    """
    local_square = np.sum(v.power(2))
    norm = cOmm.allreduce(local_square, op=MPI.SUM)
    return norm**0.5

def _gmres1_gather_and_redistribute_vector(v, vec_range, allocated):
    """
    :param scipy.sparse.csc_matrix v: The vector to be firstly collected and then distributed.
    :param allocated: A list to show how to distribute the vector?
    :return:
    """
    v = cOmm.gather(v, root=mAster_rank)
    if rAnk == mAster_rank:
        v = np.sum(v) # still a sparse matrix, but only in mAster core now
        LEN = v.shape[0]
        VList = list()
        for i in range(sIze):
            VI = spspa.csc_matrix((v[vec_range[i], 0].T.toarray()[0],
                                   vec_range[i], [0, allocated[i]]), shape=(LEN, 1))
            VList.append(VI)
    else:
        VList = None
    v = cOmm.scatter(VList, root=mAster_rank)
    return v

class _gmres1_GROUP_gather_and_redistribute_vector(FrozenOnly):
    """
    Using grouped cores to do the gathering and re-distributing.

    :param CG: Core Group information.
    :param LVR: Local Vector Range.
    :param AVR: All Vector Range.
    :param GVR: Group Vector Range.
    :param LEN: the length of the total vector, use it to put data in sparse csc vector.
    :return:
    """
    def __init__(self, CG, LVR, AVR, GVR, LEN):
        self._CG_ = CG
        self._LVR_ = LVR
        # self._AVR_ = AVR
        self._GVR_ = GVR
        # first we send data to all group leaders.
        leaders = CG['leaders']
        gtr = dict()
        for l in leaders:
            members = [l,] + CG[l]
            if len(AVR[members[0]]) == 0:
                gtr[l] = (0, 0)
            else:
                grmin = min(AVR[members[0]])
                for m in members:
                    if len(AVR[m]) > 0:
                        grmax = max(AVR[m]) + 1
                gtr[l] = (grmin, grmax)
        self._gtr_ = gtr # group total range
        if rAnk in leaders:
            self._gvl_ = 0
            for ar in GVR:
                self._gvl_ += len(ar)
        self._len_ = len(LVR)
        self._LEN_ = LEN
        self._freeze_self_()

    def __call__(self, v):
        """

        :param v:
        :return:
        """
        leaders = self._CG_['leaders']
        # gather ...
        for dest in self._gtr_:
            i0, i1= self._gtr_[dest]
            data = v[i0:i1, 0].toarray().ravel()
            if rAnk == dest:
                DATA = np.empty((sIze, self._gvl_), dtype='d')
            else:
                DATA = None
            cOmm.Gather(data, [DATA, MPI.DOUBLE], root=dest)
            if rAnk == dest:
                dt = np.sum(DATA, axis=0)
            else:
                pass
        # distribute ...
        if rAnk in leaders:
            # print(rAnk, dt, self._GVR_, self._CG_[rAnk])
            # get distributed v for self.
            IS = self._GVR_[0].start
            for k, gvr in enumerate(self._GVR_):
                i0, i1 = gvr.start, gvr.stop
                v = dt[i0-IS: i1-IS]
                if k == 0:
                    f = v
                else:
                    m = self._CG_[rAnk][k-1]
                    cOmm.Send([v, MPI.DOUBLE], dest=m, tag=rAnk)
        else:
            f = np.empty(self._len_, dtype='d')
            cOmm.Recv([f, MPI.DOUBLE], source=self._CG_['my leader'], tag=self._CG_['my leader'])
        # put f into sparse csc vector ...
        f = spspa.csc_matrix((f, self._LVR_, [0, self._len_]), shape=[self._LEN_, 1])
        # ...
        return f







class _gmres2_VIP(FrozenOnly):
    """
    Do vector inner product for gmres3 routine.
    """
    def __init__(self, AA):
        SR = AA.shared_rows
        LS_idx = list(np.argwhere(SR).ravel())

        # if rAnk == 0:
        #     print(rAnk, SR)
        # AS_idx = cOmm.gather(LS_idx, root=mAster_rank)
        indices = AA.nonempty_rows
        NS_idx = list()
        for i in indices:
            if i in LS_idx:
                pass
            else:
                NS_idx.append(i)
        assert len(set(NS_idx + LS_idx)) == len(indices), "indices division wrong."
        self._NS_idx_ = NS_idx
        self._LS_idx_ = LS_idx
        self._g_dofs_ = AA.shape[0]
        self._len_LSi_ = len(self._LS_idx_)
        self.DO_reset_cache()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._last_j_ = None
        self._last_SS_w1_cache_ = None
        self._last_NS_w1_cache_ = None
        self._Vi_cache_ = dict()

    def __call__(self, v1, v2=None, j=None, i=None):
        """
        Compute the inner product between v1 and v2.

        :param v1:
        :param v2:
        :return:
        """
        if v2 is None: # then (v1, v1)
            NS_v = v1[self._NS_idx_].data
            NS_v = np.sum(NS_v**2)

            LS_v = spspa.csc_matrix((v1[self._LS_idx_,0].T.toarray()[0],
                                     self._LS_idx_, [0, self._len_LSi_]),
                                    shape=(self._g_dofs_,1))
            LS_v = cOmm.gather(LS_v, root=sEcretary_rank)
            if rAnk == sEcretary_rank:
                LS_v = np.sum(np.sum(LS_v).data**2)
                NS_v += LS_v
            return cOmm.allreduce(NS_v, op=MPI.SUM)

        else:
            if j == self._last_j_:
                SS_w1 = self._last_SS_w1_cache_
                NS_w1 = self._last_NS_w1_cache_
            else:
                self._last_j_ = j
                SS_w1 = spspa.csc_matrix((v1[self._LS_idx_,0].T.toarray()[0],
                                         self._LS_idx_, [0, self._len_LSi_]),
                                        shape=(self._g_dofs_,1))
                SS_w1 = cOmm.gather(SS_w1, root=mAster_rank)
                if rAnk == mAster_rank: SS_w1 = np.sum(SS_w1)
                NS_w1 = v1[self._NS_idx_]

                self._last_SS_w1_cache_ = SS_w1 # only in mAster_rank, it is not None
                self._last_NS_w1_cache_ = NS_w1

            if i in self._Vi_cache_:
                SS_vi, NS_vi = self._Vi_cache_[i]
            else:
                SS_vi = spspa.csc_matrix((v2[self._LS_idx_,0].T.toarray()[0],
                                         self._LS_idx_, [0, self._len_LSi_]),
                                        shape=(self._g_dofs_,1))
                SS_vi = cOmm.gather(SS_vi, root=mAster_rank)
                if rAnk == mAster_rank: SS_vi = np.sum(SS_vi)
                NS_vi = v2[self._NS_idx_]
                self._Vi_cache_[i] = [SS_vi, NS_vi]

            NS = np.sum(NS_w1.multiply(NS_vi))
            if rAnk == mAster_rank:
                NS += np.sum(SS_w1.multiply(SS_vi))
            return cOmm.allreduce(NS, op=MPI.SUM)

            # NS = cOmm.reduce(NS, root=mAster_rank, op=MPI.SUM)
            # return cOmm.bcast(NS, root=mAster_rank)

class _gmres2_CAD(FrozenOnly):
    """Collection and distribution."""
    def __init__(self, AA):
        LNC = AA.nonempty_columns
        self._ANC_ = cOmm.gather(LNC, root=mAster_rank)
        self._total_dofs_ = AA.shape[1]
        self._LI_ = LNC
        self._nnz_ = len(self._LI_)
        self._freeze_self_()

    def __call__(self, v):
        """

        :param v: the vector to be collected and distributed.
        :return:
        """
        v = cOmm.gather(v, root=mAster_rank)
        if rAnk == mAster_rank:
            v = np.sum(v, axis=0)
            V = list()
            for i, nnc in enumerate(self._ANC_):
                if v.__class__.__name__ == 'ndarray':
                    V.append(v[nnc])
                else:
                    V.append(v[nnc, 0].T.toarray()[0])
        else:
            V = None
        v = cOmm.scatter(V, root=mAster_rank)
        v = spspa.csc_matrix((v, self._LI_, [0, self._nnz_]), shape=(self._total_dofs_, 1))
        return v

def gmres2(AA, bb, X0, restart=100, maxiter=1000, tol=1e-4):
    """

    :param AA:
    :param bb:
    :param X0:
    :param restart:
    :param maxiter:
    :param tol:
    :return:
    """
    assert AA.__class__.__name__ == 'GlobalMatrix'
    assert bb.__class__.__name__ == 'GlobalVector'
    assert X0.__class__.__name__ == 'DistributedVector'
    assert maxiter >= 1, "maxiter must be >= 1."
    assert restart >= 3, "restart must be >= 3."
    assert tol > 0, "tol must be > 0."
    bb.DO_resemble_row_distribution_of(AA) # important, after this, we can do f - A @ x0
    A = AA.M
    f = bb.V
    x0 = X0.V

    VIP = _gmres2_VIP(AA)
    CAD = _gmres2_CAD(AA)

    ITER = 0
    BETA = None
    while 1:
        r0 = f - A @ x0 # csc vector
        beta = VIP(r0) ** 0.5

        # if rAnk == 0:
        #     print('gmres2:', ITER, 'restart:', restart, 'error:', beta, flush=True)

        # check stop iteration or not ...
        if BETA is None: BETA = [beta,]
        if len(BETA) > 20: BETA = BETA[-5:]
        BETA.append(beta)
        stop_iteration, info = ___gmres_stop_criterion___(tol, ITER, maxiter, BETA)
        if stop_iteration: break
        # ...
        v0 = r0 / beta
        Vm = [v0,]
        if rAnk == sEcretary_rank: Hm = spspa.lil_matrix((restart+1, restart))

        for j in range(restart):
            vj = CAD(Vm[j])
            wj = A @ vj

            hij_vi = None
            for i in range(0, j + 1):
                Hm_ij = VIP(wj, Vm[i], j, i)
                if hij_vi is None:
                    # noinspection PyUnresolvedReferences
                    hij_vi = Hm_ij * Vm[i]
                else:
                    hij_vi += Hm_ij * Vm[i]
                if rAnk == sEcretary_rank:
                    Hm[i, j] = Hm_ij
            hat_v_jp1 = wj - hij_vi
            Hm_j1_j = VIP(hat_v_jp1)**0.5
            if rAnk == sEcretary_rank: Hm[j+1,j] = Hm_j1_j
            if j < restart-1:
                Vm.append(hat_v_jp1 / Hm_j1_j)

        VIP.DO_reset_cache() # clear cache, make it read for next

        if rAnk == sEcretary_rank:
            Hm = Hm.tocsr()
            HmT = Hm.T
            ls_A = HmT @ Hm
            ls_b = HmT[:,0] * beta
            ym = spspalinalg.spsolve(ls_A, ls_b)
            del HmT, ls_A, ls_b
        else:
            ym = np.empty(restart, dtype='d')
        cOmm.Bcast([ym, MPI.DOUBLE], root=sEcretary_rank)
        Vm = spspa.hstack(Vm)
        x0 += CAD(Vm @ ym)
        ITER += 1

    if info < 0:
        raise LinerSystemSolverDivergenceError(
            f"gmres2 diverges after {ITER} iterations with error reaching {beta}.")

    return x0, info, beta, ITER
