# -*- coding: utf-8 -*-
from components.exceptions import LocalCochainShapeError
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_ColumnVector
from scipy.sparse import lil_matrix, csr_matrix, csc_matrix
from tools.miLinearAlgebra.dataStructures.vectors.distributed.main import DistributedVector
from components.freeze.base import FrozenOnly
from root.config.main import np, RANK, MASTER_RANK, COMM
from objects.CSCG.base.forms.standard.cochain.dofwise import CSCG_SF_Cochain_DofWise





class CSCG_Standard_Form_Cochain_BASE(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._local_ = None
        self._dofwise_ = None
        self._freeze_self_()

    @property
    def dofwise(self):
        """f.cochain.dofwise[i] will return the cochain for the dof #i in all cores."""
        if self._dofwise_ is None:
            self._dofwise_ = CSCG_SF_Cochain_DofWise(self)
        return self._dofwise_

    @property
    def array(self):
        """return the local cochain as a 2d array."""
        if len(self) == 0:
            return None
        else:
            return np.array(list(self.local.values()))

    @property
    def EWC(self):
        """Return the cochain as an Element-Wise-Cached vector.

        Notice that if we have changed the local cochain, the EWC will also change because we make the vector in real
        time.

        """
        ewc = EWC_ColumnVector(self._sf_.mesh.elements, self.___PRIVATE_local_call___)
        ewc.gathering_matrix = self._sf_
        return ewc


    def ___PRIVATE_local_call___(self, i):
        return csr_matrix(self.local[i]).T


    def ___DO_gather_local_to_master___(self):
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
            rn, loc_ind = self._sf_.mesh.do.find.region_name_and_local_indices_of_element(i)

            RN_LI_dict[i] = rn + '=|=' + str(loc_ind)

        RN_LI_dict = COMM.gather(RN_LI_dict, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            RID = dict()
            for rid in RN_LI_dict:
                RID.update(rid)
        del RN_LI_dict

        LOCAL = self.___DO_gather_local_to_master___()
        if RANK == MASTER_RANK:

            RW_LOCAL = dict()

            for i in range(self._sf_.mesh.elements.GLOBAL_num):
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

        for i in self._sf_.mesh.elements:
            rn, loc_ind = self._sf_.mesh.do.find.region_name_and_local_indices_of_element(i)
            dict_key = rn + '=|=' + str(loc_ind)
            LOC_COCHAIN[i] = RW_LI_COCHAIN[dict_key]

        self.local = LOC_COCHAIN





    @property
    def globe(self):
        GM = self._sf_.numbering.gathering
        globe = lil_matrix((1, self._sf_.num.global_dofs))
        for i in GM: # go through all local elements
            globe[0, GM[i].full_vector] = self.local[i]
        globe = globe.tocsr().T

        if self._sf_.whether.hybrid or (self._sf_.ndim == self._sf_.k):
            # clearly, for volume form or hybrid form, we can make a distributed vector directly
            return DistributedVector(globe)

        else:
            # otherwise, we make a master dominating distributed vector.
            GLOBE = COMM.gather(globe, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                measure = np.zeros(self._sf_.num.global_dofs, dtype=int)
                for G in GLOBE:
                    indices = G.indices
                    measure[indices] += 1

                measure[measure==0] = 1
                # noinspection PyUnresolvedReferences
                _____ = np.sum(GLOBE).toarray().ravel() / measure
                globe = csr_matrix(_____).T

            else:
                globe = csc_matrix((self._sf_.num.global_dofs, 1))

            GDV = DistributedVector(globe)
            assert GDV.whether.master_dominating
            return GDV

    @globe.setter
    def globe(self, globe):
        """
        This process is complex, but it makes sure that the distribution is correct for all cases.

        :param globe:
        :return:
        """
        if globe.__class__.__name__ == 'DistributedVector':
            assert globe.V.shape == (self._sf_.num.global_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if globe.whether.master_dominating:
                # no need to gather
                VV = globe.V.T.toarray()[0]
            else:
                V = globe.V
                V = COMM.gather(V, root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    VV = np.empty((self._sf_.num.global_dofs,))
                    for v in V:
                        indices = v.indices
                        data = v.data
                        VV[indices] = data
            # distribute vector to individual cores ...
            local_range = self._sf_.numbering.gathering.local_range
            local_range = COMM.gather(local_range, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                TO_BE_SENT = list()
                for lr in local_range:
                    if lr == tuple():
                        to_be_sent = None
                    else:
                        # noinspection PyUnboundLocalVariable
                        to_be_sent = csc_matrix(
                            (VV[lr[0]:lr[1]], range(lr[0],lr[1]), [0, lr[1]-lr[0]]),
                            shape=(self._sf_.num.global_dofs, 1))
                    TO_BE_SENT.append(to_be_sent)
            else:
                TO_BE_SENT = None
            TO_BE_SENT = COMM.scatter(TO_BE_SENT, root=MASTER_RANK)
            # distribute to local cochain ...
            local = dict()
            GM = self._sf_.numbering.gathering
            for i in GM: # go through all local elements
                idx = GM[i].full_vector
                local[i] = TO_BE_SENT[idx].toarray().ravel()
            self.local = local

        elif globe.__class__.__name__ == 'LocallyFullVector':
            V = globe.V # V already be 1-d array.
            local = dict()
            GM = self._sf_.numbering.gathering
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


    # ------------- DEPENDENT PROPERTIES (MAJOR): When set, clear BRANCHES by set _branches_ to None -------------------
    @property
    def local(self):
        """The local cochain. Must be full. So all local mesh elements must have their local cochains!

        :return: A dict whose keys are local element indices and values are cochains (1-d arrays) in corresponding
            elements.
        :rtype: Dict[int, numpy.ndarray]
        """
        return self._local_

    @local.setter
    def local(self, local):
        numOfElements = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.num.basis
        try:
            assert isinstance(local, dict), \
                f"local cochain needs to be a dict, now it is a {local.__class__.__name__}."
            assert len(local) == numOfElements, \
                "local cochain has to contain cochains for all local mesh elements."
            for i, j in zip(self._sf_.mesh.elements.indices, local):
                assert np.shape(local[i]) == (numOfBasis,), \
                    f"local[{i}] shape = {np.shape(local[i])} wrong. " \
                    f"It needs to be {(numOfBasis,)}."
                assert i == j, f"mesh element index sequence is wrong."

        except AssertionError:
            raise LocalCochainShapeError("Cannot set local cochain.")

        self._local_ = local


    #--------------- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set, update local ------


    #=================================== ABOVE =========================================================================

