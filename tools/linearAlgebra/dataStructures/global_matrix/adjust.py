# -*- coding: utf-8 -*-


from screws.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK



class ___GM_ADJUST___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    # PRIVATE, very low efficiency, please do not use --------------------- BELOW ------------------
    def ___PRIVATE_clear_row___(self, r):
        """
        Make row #i all zero.

        This is done in all cores such that the row is for sure zero then.

        :param r:
        :return: None
        """
        if self._gm_.mtype != 'lil':
            self._gm_._M_ = self._gm_._M_.tolil() # this will
        self._gm_._M_[r, :] = 0

    def ___PRIVATE_set_value___(self, i, j, value):
        """
        Set M[i,j] = v only in master core, set M[i,j] = 0 in all other cores.

        :param i:
        :param j:
        :param value:
        :return:
        """
        if self._gm_.mtype != 'lil':
            self._gm_._M_ = self._gm_._M_.tolil()

        if RANK == MASTER_RANK:
            self._gm_._M_[i, j] = value
        else:
            self._gm_._M_[i, j] = 0
    # -------------------------- ABOVE -------------------------------------------------------------