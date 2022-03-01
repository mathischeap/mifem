# -*- coding: utf-8 -*-
"""
"""
from root.config import *
from screws.frozen import FrozenOnly
from screws.exceptions import LocalCochainShapeError
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector
from scipy.sparse import lil_matrix, csr_matrix, csc_matrix
from tools.linear_algebra.data_structures.global_matrix.main import DistributedVector
from scipy import sparse as spspa


# noinspection PyUnresolvedReferences
class CSCG_Standard_Form:
    """"""
    @property
    def NUM_basis(self):
        """(int) Return a int which represent the number of basis function one element has."""
        return self._NUM_basis_

    @property
    def NUM_basis_components(self):
        """
        (Tuple[int]) Return a tuple of integers. If it is 0- or 3-form, then it is like (x,),
        and when it is 1- or 2-form, then it is like (x,y,z) because it represents a vector. But after all,
        sum(NUM_basis_components) == NUM_basis.
        """
        return self._NUM_basis_components_

    @property
    def NUM_dofs(self):
        """(int) Return number of dofs in this core."""
        return self.numbering.num_of_dofs_in_this_core

    @property
    def GLOBAL_num_dofs(self):
        """(int) Return number of total dofs."""
        return self.numbering.gathering.GLOBAL_num_dofs


    @property
    def orientation(self):
        return self._orientation_


    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def error(self):
        return self._error_

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def operators(self):
        return self._operators_

    @property
    def visualize(self):
        return self._visualize_



    @property
    def IS_hybrid(self):
        return self._IS_hybrid_

    @property
    def IS_inner_oriented(self):
        return True if self._orientation_ == 'inner' else False


    @property
    def IS_volume_form(self):
        return self.ndim == self.k

    @property
    def do(self):
        """If it has too many do methods, we group them in to do."""
        return self._DO_








class CSCG_Standard_Form_Cochain_BASE(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._local_ = None # this is a key property, should not reset it.
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        pass


    @property
    def EWC(self):
        """Return the cochain as an Element-Wise-Cached vector.

        Notice that if we have changed the local cochain, the EWC will also change because we make the vector in real
        time.

        """
        ewc = EWC_ColumnVector(self._sf_.mesh.elements, self.___PRIVATE_local_call___)
        ewc.gathering_matrix = self._sf_
        return ewc


    def ___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___(self):
        """make it regions-wise-element-local-indexed, thus we can save it and when read a form, we can always have the
        correct local cochain allocated even element numbering is different.
        """
        assert self.local is not None, "I have no local cochain!"

        RN_LI_dict = dict()
        for i in self.local:
            rn, loc_ind = self._sf_.mesh.do.find.region_name_and_local_indices_of_element(i)

            RN_LI_dict[i] = rn + '=|=' + str(loc_ind)

        RN_LI_dict = cOmm.gather(RN_LI_dict, root=mAster_rank)
        if rAnk == mAster_rank:
            RID = dict()
            for rid in RN_LI_dict:
                RID.update(rid)
        del RN_LI_dict

        LOCAL = self.DO_gather_local_to_master()
        if rAnk == mAster_rank:

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


    def ___PRIVATE_local_call___(self, i):
        return csr_matrix(self.local[i]).T

    @property
    def globe(self):
        GM = self._sf_.numbering.gathering
        globe = lil_matrix((1, self._sf_.GLOBAL_num_dofs))
        for i in GM: # go through all local elements
            globe[0, GM[i].full_vector] = self.local[i]
        globe = globe.tocsr().T

        if self._sf_.IS_hybrid or (self._sf_.ndim == self._sf_.k):
            # clearly, for volume form or hybrid form, we can make a distributed vector directly
            return DistributedVector(globe)

        else:
            # otherwise, we make a master dominating distributed vector.
            GLOBE = cOmm.gather(globe, root=mAster_rank)

            if rAnk == mAster_rank:
                measure = np.zeros(self._sf_.GLOBAL_num_dofs, dtype=int)
                for G in GLOBE:
                    indices = G.indices
                    measure[indices] += 1

                measure[measure==0] = 1
                # noinspection PyUnresolvedReferences
                _____ = np.sum(GLOBE).toarray().ravel() / measure
                globe = csr_matrix(_____).T

            else:
                globe = csc_matrix((self._sf_.GLOBAL_num_dofs,1))

            GDV = DistributedVector(globe)
            assert GDV.IS_master_dominating
            return GDV

    @globe.setter
    def globe(self, globe):
        """
        This process is complex, but it makes sure that the distribution is correct for all cases.

        :param globe:
        :return:
        """
        if globe.__class__.__name__ == 'DistributedVector':
            assert globe.V.shape == (self._sf_.GLOBAL_num_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if globe.IS_master_dominating:
                # no need to gather
                VV = globe.V.T.toarray()[0]
            else:
                V = globe.V
                V = cOmm.gather(V, root=mAster_rank)
                if rAnk == mAster_rank:
                    VV = np.empty((self._sf_.GLOBAL_num_dofs,))
                    for v in V:
                        indices = v.indices
                        data = v.data
                        VV[indices] = data
            # distribute vector to individual cores ...
            local_range = self._sf_.numbering.gathering.local_range
            local_range = cOmm.gather(local_range, root=mAster_rank)
            if rAnk == mAster_rank:
                TO_BE_SENT = list()
                for lr in local_range:
                    if lr == tuple():
                        to_be_sent = None
                    else:
                        # noinspection PyUnboundLocalVariable
                        to_be_sent = spspa.csc_matrix(
                            (VV[lr[0]:lr[1]], range(lr[0],lr[1]), [0, lr[1]-lr[0]]),
                            shape=(self._sf_.GLOBAL_num_dofs, 1))
                    TO_BE_SENT.append(to_be_sent)
            else:
                TO_BE_SENT = None
            TO_BE_SENT = cOmm.scatter(TO_BE_SENT, root=mAster_rank)
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


    def DO_gather_local_to_master(self):
        """Do what the method name says."""
        local = cOmm.gather(self.local, root=mAster_rank)
        if rAnk == mAster_rank:
            LOCAL = dict()
            for li in local:
                if li is not None:
                    LOCAL.update(li)
            return LOCAL


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
        """The local cochain.

        :return: A dict whose keys are local element indices and values are cochains (1-d arrays) in corresponding
            elements.
        :rtype: Dict[int, numpy.ndarray]
        """
        return self._local_

    @local.setter
    def local(self, local):
        numOfElements = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.NUM_basis
        try:
            assert isinstance(local, dict), f"local cochain needs to be a dict, now it is a {local.__class__.__name__}."
            assert len(local) == numOfElements, "local cochain has to contain cochain for all local elements."
            for i in self._sf_.mesh.elements.indices:
                assert np.shape(local[i]) == (numOfBasis,), f"local[{i}] shape = {np.shape(local[i])} wrong. " \
                                                            f"It needs to be {(numOfBasis,)}."
        except AssertionError:
            raise LocalCochainShapeError("Cannot set local cochain.")
        self._local_ = local


    #--------------- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set, update local ------


    #=================================== ABOVE =========================================================================






class CSCG_Standard_Form_Coboundary_BASE(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._next_form_ = None
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._incidenceMatrix_ = None

    @property
    def incidence_matrix(self):
        raise NotImplementedError()

    def ___PRIVATE_next_class___(self):
        raise NotImplementedError()

    def __call__(self):
        """
        When we call the coboundary object, we do the ``coboundary`` process; let ``self`` be a ``k``-form,
        it returns a ``(k+1)``-form.

        :return: A new standard ``(k+1)``-form.
        :raise AssertionError: If ``self`` has no cochain.
        """
        assert self._sf_.cochain.local is not None, "I need a cochain to perform coboundary."
        nextFmClass = self.___PRIVATE_next_class___()
        nextFmInstance = nextFmClass(
            self._sf_.mesh, self._sf_.space,
            is_hybrid = self._sf_.IS_hybrid,
            numbering_parameters = self._sf_.numbering._numbering_parameters_,
            name = 'd(' + self._sf_.standard_properties.name + ')'
        )
        selfCochain = self._sf_.cochain.local
        nextCochain = dict()
        incidence_matrix = self.incidence_matrix
        for i in self._sf_.mesh.elements:
            nextCochain[i] = incidence_matrix[i] @ selfCochain[i]
        nextFmInstance.cochain.local = nextCochain
        return nextFmInstance

    @property
    def cochain_local(self):
        """
        Return the local cochain (not the form) of its coboundary.

        :return:
        """
        selfCochain = self._sf_.cochain.local
        nextCochain = dict()
        incidence_matrix = self.incidence_matrix
        for i in self._sf_.mesh.elements:
            nextCochain[i] = incidence_matrix[i] @ selfCochain[i]
        return nextCochain
