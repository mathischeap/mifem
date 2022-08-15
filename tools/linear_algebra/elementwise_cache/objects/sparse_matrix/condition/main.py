# -*- coding: utf-8 -*-


from screws.freeze.base import FrozenOnly


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