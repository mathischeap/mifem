# -*- coding: utf-8 -*-
import numpy as np

from root.config.main import COMM, MASTER_RANK, RANK
from components.freeze.base import FrozenOnly


class EWC_SpaMat_Condition(FrozenOnly):
    """"""
    def __init__(self, spa_mat):
        """"""
        self._spa_mat_ = spa_mat
        self._freeze_self_()

    @property
    def condition_number(self):
        """{float} : the condition number of the assembled matrix."""
        A = self._spa_mat_.assembled
        condition_number = A.condition.condition_number
        return condition_number

    @property
    def pseudo_sparsity(self):
        """We will not assemble the matrix. The sparsity is equal to the mean sparsity of all
        elements. Which of course is not correct. But it at least tells something. That is why
        it is called `pseudo_sparsity`.

        """

        local_sparsity = list()
        for i in self._spa_mat_:
            spa_mat = self._spa_mat_[i]

            local_sparsity.append(1 - spa_mat.nnz / np.prod(spa_mat.shape))

        if len(local_sparsity) > 0:
            local_sparsity = np.sum(local_sparsity) / len(local_sparsity)
        else:
            local_sparsity = None

        local_sparsity = COMM.gather(local_sparsity, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            valid = list()
            for ls in local_sparsity:
                if ls is not None:
                    valid.append(ls)

            local_sparsity = sum(valid) / len(valid)

        return COMM.bcast(local_sparsity, root=MASTER_RANK)