# -*- coding: utf-8 -*-
from root.config.main import np, RANK, MASTER_RANK, COMM, SAFE_MODE
from tools.linearAlgebra.dataStructures.global_matrix.main import DistributedVector
from scipy import sparse as spspa
from screws.exceptions import LocalCochainShapeError
from scipy.sparse import lil_matrix, csr_matrix, csc_matrix
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_ColumnVector
from screws.freeze.base import FrozenOnly


class CSCG_Trace_Form_Cochain_BASE(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._local_ = None
        self._local_TEW_ = None
        self._freeze_self_()

    def RESET_cache(self):
        self._local_TEW_ = None

    @property
    def EWC(self):
        """Return the cochain as an Element-Wise-Cached vector.

        Notice that if we have changed the local cochain, the EWC will also change because we make the vector in real
        time.

        """
        ewc = EWC_ColumnVector(self._tf_.mesh.elements, self.___PRIVATE_local_call___)
        ewc.gathering_matrix = self._tf_
        return ewc

    def ___PRIVATE_local_call___(self, i):
        return csr_matrix(self.local[i]).T

    def ___PRIVATE_gather_local_to_master___(self):
        """Do what the method name says."""
        local = COMM.gather(self.local, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            LOCAL = dict()
            for li in local:
                if li is not None:
                    LOCAL.update(li)
            return LOCAL

    def ___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___(self):
        """make it regions-wise-element-local-indexed, thus we can save it and when read a form, we can always have the
        correct local cochain allocated even element numbering is different.
        """
        assert self.local is not None, "I have no local cochain!"

        RN_LI_dict = dict()
        for i in self.local:
            rn, loc_ind = self._tf_.mesh.do.find.region_name_and_local_indices_of_element(i)

            RN_LI_dict[i] = rn + '=|=' + str(loc_ind)

        RN_LI_dict = COMM.gather(RN_LI_dict, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            RID = dict()
            for rid in RN_LI_dict:
                RID.update(rid)
        del RN_LI_dict

        LOCAL = self.___PRIVATE_gather_local_to_master___()
        if RANK == MASTER_RANK:

            RW_LOCAL = dict()

            for i in range(self._tf_.mesh.elements.GLOBAL_num):
                assert i in LOCAL, "something is wrong."
                # noinspection PyUnboundLocalVariable
                assert i in RID, "something is wrong."

                rn_loc_ind = RID[i]
                RW_LOCAL[rn_loc_ind] = LOCAL[i]

            return RW_LOCAL

    def ___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(self, RW_LI_COCHAIN):
        """When we have the Region-wised local index grouped cochain, we can use this method to distribute it to local
        cochain. The regions-wise grouped cochain must be a full cochain in all cores.

        :param RW_LI_COCHAIN:
        :return:
        """

        LOC_COCHAIN = dict()

        for i in self._tf_.mesh.elements:
            rn, loc_ind = self._tf_.mesh.do.find.region_name_and_local_indices_of_element(i)
            dict_key = rn + '=|=' + str(loc_ind)
            LOC_COCHAIN[i] = RW_LI_COCHAIN[dict_key]

        self.local = LOC_COCHAIN

    @property
    def globe(self):
        """Global cochain. As trace elements almost are always shared by cores, so we
        cannot make it a general`DistributedVector`; we can only make it a master dominating one.
        """
        GM = self._tf_.numbering.gathering
        globe = lil_matrix((1, self._tf_.GLOBAL_num_dofs))
        for i in GM: # go through all local elements
            globe[0, GM[i].full_vector] = self.local[i]

        globe = globe.tocsr().T
        GLOBE = COMM.gather(globe, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            measure = np.zeros(self._tf_.GLOBAL_num_dofs, dtype=int)
            for G in GLOBE:
                indices = G.indices
                measure[indices] += 1

            measure[measure == 0] = 1
            # noinspection PyUnresolvedReferences
            _____ = np.sum(GLOBE).toarray().ravel() / measure
            globe = csr_matrix(_____).T

        else:
            globe = csc_matrix((self._tf_.GLOBAL_num_dofs, 1))

        GDV = DistributedVector(globe)
        assert GDV.IS.master_dominating
        return GDV

    @globe.setter
    def globe(self, globe):
        if globe.__class__.__name__ == 'DistributedVector':
            assert globe.V.shape == (self._tf_.GLOBAL_num_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if globe.IS_master_dominating:
                # no need to gather
                VV = globe.V.T.toarray()[0]
            else:
                V = globe.V
                V = COMM.gather(V, root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    VV = np.empty((self._tf_.GLOBAL_num_dofs,))
                    for v in V:
                        indices = v.indices
                        data = v.data
                        VV[indices] = data
            # distribute vector to individual cores ...
            local_range = self._tf_.numbering.gathering.local_range
            local_range = COMM.gather(local_range, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                TO_BE_SENT = list()
                for lr in local_range:
                    if lr == tuple():
                        to_be_sent = None
                    else:
                        # noinspection PyUnboundLocalVariable
                        to_be_sent = spspa.csc_matrix(
                            (VV[lr[0]:lr[1]], range(lr[0],lr[1]), [0, lr[1]-lr[0]]),
                            shape=(self._tf_.GLOBAL_num_dofs, 1))
                    TO_BE_SENT.append(to_be_sent)
            else:
                TO_BE_SENT = None
            TO_BE_SENT = COMM.scatter(TO_BE_SENT, root=MASTER_RANK)
            # distribute to local cochain ...
            local = dict()
            GM = self._tf_.numbering.gathering
            for i in GM: # go through all local elements
                idx = GM[i].full_vector
                local[i] = TO_BE_SENT[idx].toarray().ravel()
            self.local = local

        elif globe.__class__.__name__ == 'LocallyFullVector':
            V = globe.V # V already be 1-d array.
            local = dict()
            GM = self._tf_.numbering.gathering
            for i in GM:  # go through all local elements
                idx = GM[i].full_vector
                local[i] = V[idx]
            self.local = local

        else:
            raise Exception(f"Can not set cochain from {globe}.")



    def __getitem__(self, item):
        return self.local[item]

    def __contains__(self, item):
        return item in self.local

    def __iter__(self):
        for i in self.local:
            yield i

    def __len__(self):
        return len(self.local)


    #--DEPENDENT PROPERTIES (MAJOR): When set local, clear BRANCHES by set branches to None ------
    @property
    def local(self):
        """
        The local cochain. Must be full (all local mesh elements have their local cochains.)

        :return: A dict whose keys are local element indices and values are cochain in corresponding elements.
        :rtype: Dict[int, numpy.ndarray]
        """
        return self._local_

    @local.setter
    def local(self, local):
        numOfElements = self._tf_.mesh.elements.num
        numOfBasis = self._tf_.num.basis
        try:
            assert isinstance(local, dict)
            assert len(local) == numOfElements
            for i in self._tf_.mesh.elements:
                assert np.shape(local[i]) == (numOfBasis,)
        except AssertionError:
            raise LocalCochainShapeError()

        self.RESET_cache()
        self._local_ = local

    #--DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set below, update local ------
    @property
    def local_TEW(self):
        """
        The local cochain in trace element.

        TEW stands for Trace-Element-Wise.

        :return: A dict whose keys are trace element names and values are cochain in corresponding trace elements.
        :rtype: Dict[str, numpy.ndarray]
        """
        # this is important: do not use ``local_TEW`` or ``local``.
        if self._local_TEW_ is None and self._local_ is not None:
            self.___local_2_local_TEW___()
        return self._local_TEW_

    @local_TEW.setter
    def local_TEW(self, local_TEW):
        """
        :param local_TEW:
        :return:
        """
        numOfBasis = self._tf_.num.basis_onside
        try:
            assert isinstance(local_TEW, dict)
            for key in local_TEW:
                te = self._tf_.mesh.trace.elements[key]
                rs = te.CHARACTERISTIC_side
                assert local_TEW[key].shape == (numOfBasis[rs],)
        except AssertionError:
            raise LocalCochainShapeError()

        self.RESET_cache()
        self._local_TEW_ = local_TEW
        self.___local_TEW_2_local___()

    def ___local_TEW_2_local___(self):
        """"""
        MAP = self._tf_.mesh.trace.elements.map
        if SAFE_MODE:
            for i in MAP:
                for key in MAP[i]:
                    assert key in self._local_TEW_, "'local_TEW' is not full."
        local = dict()
        for i in MAP:
            local[i] = list()
            for key in MAP[i]:
                local[i].append(self._local_TEW_[key])
            local[i] = np.concatenate(local[i])
            assert local[i].shape == (self._tf_.num.basis,)
        self._local_ = local

    def ___local_2_local_TEW___(self):
        return NotImplementedError()

    #=================================== ABOVE -> branch 1 =======================================