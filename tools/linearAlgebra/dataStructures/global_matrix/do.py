# -*- coding: utf-8 -*-


from tools.linearAlgebra.dataStructures.global_matrix.helpers.split_A import ___split_A___

from screws.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK, COMM, TREE
import numpy as np
from scipy import sparse as spspa

class ___GM_DO___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    def claim_distribution_pattern(self):
        """
        We parse the structure of M and classify the global matrix into 'row', 'column' or False.

        :return:
        """
        IS_row_major = self._gm_.___PRIVATE_check_row_major___()
        if IS_row_major:
            self._gm_.IS.regularly_distributed = 'row'
            return

        IS_column_major = self._gm_.___PRIVATE_check_col_major___()
        if IS_column_major:
            self._gm_.IS.regularly_distributed = 'column'
            return

        self._gm_.IS.regularly_distributed = False


    def gather_M_to_core(self, core=None,
        clean_local=False, splitting_factor=50000000):
        """
        We gather M to one core, leave M to be of all zero elements in all other cores.

        Unlike the gathering method for vector, this gathering method will damage the original data
        if ``clean_local`` is True.

        Do not easily use this, this can be very slow and can damage the data or cause
        failure for MPI.

        :param core: The core that recv all matrix elements.
        :param bool clean_local:
        :param int splitting_factor: When the matrix is too big, the MPI may fail (can not pass
            object larger than 2GB). So we need to split it. The criterion is: If the number of
            non-zero values are larger than ``splitting_factor`` we do the splitting. Remember,
            splitting will only be used when we do ``clean_local`` after gathering.
        :return:
        """
        SELF  = self._gm_
        if core is None: core = MASTER_RANK

        if core == MASTER_RANK and SELF.IS.master_dominating:
            # In this case, we directly return.
            return SELF.M

        splitting_factor = int(splitting_factor)
        assert splitting_factor > 0, f"splitting_factor={splitting_factor} wrong, must be > 0."

        if clean_local:
            if SELF.GLOBAL_approx_nnz < splitting_factor:
                M = COMM.gather(SELF.M, root=core)
                if RANK == core:
                    SELF._M_ = np.sum(M)
                else:
                    SELF._M_ = spspa.csr_matrix(SELF.shape)
            else:
                assert core == MASTER_RANK, "This routine only work for root=master yet!"
                tree = TREE(2)
                for Hi in tree:  # combine `A` to master core, `A` will be cleaned.
                    if Hi is None:
                        pass

                    elif Hi[0] == 'send':
                        if SELF.nnz <= splitting_factor:
                            ST, SL = 1, None
                        else:
                            if SELF.mtype == 'csr':
                                PS = SELF.shape[0]
                            elif SELF.mtype == 'csc':
                                PS = SELF.shape[1]
                            else:
                                raise Exception()
                            ST, SL = ___split_A___(SELF.M.indptr, splitting_factor, PS)

                        COMM.send(ST, **Hi[1])
                        if ST == 1:
                            COMM.send(SELF.M, **Hi[1])
                            SELF._M_ = spspa.csr_matrix(SELF.shape)
                        else:
                            for t in range(ST):
                                if SELF.mtype == 'csr':
                                    tbs = SELF.M[SL[t]:SL[t+1], :]
                                elif SELF.mtype == 'csc':
                                    tbs = SELF.M[:, SL[t]:SL[t+1]]
                                else:
                                    raise Exception()
                                COMM.send(tbs, **Hi[1])
                            SELF._M_ = spspa.csr_matrix(SELF.shape)

                    elif Hi[0] == 'recv':
                        RT = COMM.recv(**Hi[1])
                        if RT == 1:
                            SELF._M_ += COMM.recv(**Hi[1])
                        else:
                            RECV = list()
                            for t in range(RT):
                                RECV.append(COMM.recv(**Hi[1]))
                            if SELF.mtype == 'csr':
                                SELF._M_ += spspa.vstack(RECV, format='csr')
                            elif SELF.mtype == 'csc':
                                SELF._M_ += spspa.hstack(RECV, format='csc')
                            else:
                                raise Exception()

                    else:
                        raise Exception()

            if RANK == MASTER_RANK: SELF._M_.sum_duplicates()

            return SELF._M_

        else:
            M = COMM.gather(SELF.M, root=core)
            if RANK == core:
                return np.sum(M)
            else:
                return spspa.csc_matrix(SELF.shape)