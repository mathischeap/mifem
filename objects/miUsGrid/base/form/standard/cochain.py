# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import numpy as np
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_ColumnVector
from scipy.sparse import csr_matrix
from root.config.main import RANK, MASTER_RANK, COMM
from scipy.sparse import lil_matrix, csc_matrix
from tools.miLinearAlgebra.dataStructures.vectors.distributed.main import DistributedVector


class miUsGrid_SF_CochainBase(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._local_ = None
        self._freeze_self_()

    def RESET_cache(self):
        """"""

    def __getitem__(self, item):
        return self.local[item]

    def __contains__(self, item):
        return item in self.local

    def __iter__(self):
        for i in self.local:
            yield i

    def __len__(self):
        return len(self.local)



    @property
    def globe(self):
        """"""
        GM = self._sf_.numbering.gathering
        globe = lil_matrix((1, self._sf_.num.global_dofs))
        for i in GM: # go through all local elements
            globe[0, GM[i].full_vector] = self.local[i]
        globe = globe.tocsr().T

        if self._sf_.ndim == self._sf_.k:
            # clearly, for a volume form, we can make a distributed vector directly
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
    def globe(self, glb):
        """
        This process is complex, but it makes sure that the distribution is correct for all cases.

        :param glb:
        :return:
        """
        if glb.__class__.__name__ == 'DistributedVector':
            assert glb.V.shape == (self._sf_.num.global_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if glb.whether.master_dominating:
                # no need to gather
                VV = glb.V.T.toarray()[0]
            else:
                V = glb.V
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

        elif glb.__class__.__name__ == 'LocallyFullVector':
            V = glb.V # V already be 1-d array.
            local = dict()
            GM = self._sf_.numbering.gathering
            for i in GM:  # go through all local elements
                idx = GM[i].full_vector
                local[i] = V[idx]
            self.local = local

        else:
            raise Exception(f"Can not set cochain from {glb}.")

    @property
    def EWC(self):
        """Return the cochain as an Element-Wise-Cached vector.

        Notice that if we have changed the local cochain, the EWC will also change because we make the vector in real
        time.

        """
        ewc = EWC_ColumnVector(self._sf_.mesh.elements,
                               self.___PRIVATE_local_call___,
                               cache_key_generator='no_cache')
        ewc.gathering_matrix = self._sf_
        return ewc

    def ___PRIVATE_local_call___(self, i):
        return csr_matrix(self.local[i]).T

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
        numOfElements = self._sf_.mesh.elements.num.cells
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
            raise Exception("Cannot set local cochain.")

        self.RESET_cache()
        self._local_ = local


    #--- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set, update local ------


    #======================== ABOVE =========================================================================

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
