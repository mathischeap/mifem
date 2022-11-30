# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_ColumnVector
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix
from root.config.main import RANK, MASTER_RANK, COMM
import numpy as np
from tools.miLinearAlgebra.dataStructures.vectors.distributed.main import DistributedVector

from components.exceptions import LocalCochainShapeError


class CSCG_LocalTrace_CochainBase(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._local_ = None
        self._freeze_self_()

    def RESET_cache(self):
        """"""
        pass

    @property
    def EWC(self):
        """Return the cochain as an Element-Wise-Cached vector.

        Notice that if we have changed the local cochain, the EWC will also change because we make the vector in real
        time.

        """
        ewc = EWC_ColumnVector(self._ltf_.mesh.elements, self.___PRIVATE_local_call___)
        ewc.gathering_matrix = self._ltf_.numbering.gathering
        return ewc


    def ___PRIVATE_local_call___(self, i):
        return csr_matrix(self.local[i]).T






    @property
    def globe(self):
        GM = self._ltf_.numbering.gathering
        globe = lil_matrix((1, self._ltf_.num.global_dofs))
        for i in GM: # go through all local elements
            globe[0, GM[i].full_vector] = self.local[i]
        globe = globe.tocsr().T

        # local trace is actually automatically hybrid. So it must can form a distributed vector.
        return DistributedVector(globe)


    @globe.setter
    def globe(self, globe):
        """
        This process is complex, but it makes sure that the distribution is correct for all cases.

        :param globe:
        :return:
        """
        if globe.__class__.__name__ == 'DistributedVector':
            assert globe.V.shape == (self._ltf_.num.global_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if globe.whether.master_dominating:
                # no need to gather
                VV = globe.V.T.toarray()[0]
            else:
                V = globe.V
                V = COMM.gather(V, root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    VV = np.empty((self._ltf_.num.global_dofs,))
                    for v in V:
                        indices = v.indices
                        data = v.data
                        VV[indices] = data
            # distribute vector to individual cores ...
            local_range = self._ltf_.numbering.gathering.local_range
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
                            shape=(self._ltf_.num.global_dofs, 1))
                    TO_BE_SENT.append(to_be_sent)
            else:
                TO_BE_SENT = None
            TO_BE_SENT = COMM.scatter(TO_BE_SENT, root=MASTER_RANK)
            # distribute to local cochain ...
            local = dict()
            GM = self._ltf_.numbering.gathering
            for i in GM: # go through all local elements
                idx = GM[i].full_vector
                local[i] = TO_BE_SENT[idx].toarray().ravel()
            self.local = local

        elif globe.__class__.__name__ == 'LocallyFullVector':
            V = globe.V # V already be 1-d array.
            local = dict()
            GM = self._ltf_.numbering.gathering
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
        numOfElements = self._ltf_.mesh.elements.num
        numOfBasis = self._ltf_.num.basis
        try:
            assert isinstance(local, dict), \
                f"local cochain needs to be a dict, now it is a {local.__class__.__name__}."
            assert len(local) == numOfElements, \
                "local cochain has to contain cochains for all local mesh elements."
            for i, j in zip(self._ltf_.mesh.elements.indices, local):
                assert np.shape(local[i]) == (numOfBasis,), \
                    f"local[{i}] shape = {np.shape(local[i])} wrong. " \
                    f"It needs to be {(numOfBasis,)}."
                assert i == j, f"mesh element index sequence is wrong."

        except AssertionError:
            raise LocalCochainShapeError("Cannot set local cochain.")

        self.RESET_cache()
        self._local_ = local


    #--------------- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set, update local ------


    #=================================== ABOVE =========================================================================


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
