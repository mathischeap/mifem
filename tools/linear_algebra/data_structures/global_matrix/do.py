# -*- coding: utf-8 -*-



from screws.freeze.main import FrozenOnly

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


    def gather_M_to_core(self, **kwargs):
        """A wrapper of `___PRIVATE_gather_M_to_core___` method."""
        return self._gm_.___PRIVATE_gather_M_to_core___(**kwargs)