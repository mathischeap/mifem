# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
import numpy as np
from tools.linear_algebra.gathering.vector import Gathering_Vector
from tools.linear_algebra.gathering.regular.matrix.do.find import Gathering_Matrix_FIND


class Gathering_Matrix_DO(FrozenOnly):
    """"""

    def __init__(self, GM):
        self._GM_ = GM
        self._find_ = None
        self._freeze_self_()


    def hstack(self, *args):
        """`hstack` other gathering matrices to the right of self.

        For example, we can do GM1.DO_hstack(GM2, GM3).

        :param args:
        :return:
        """
        SELF = self._GM_
        for ogm in args:
            assert ogm.__class__.__name__ == 'Gathering_Matrix', \
                f"I am stacking a {ogm.__class__.__name__}, wrong!"
            assert len(SELF) == len(ogm), \
                f"Length dis-match."

        for i in SELF:
            for ogm in args:
                assert i in ogm, "elements dis-match."

        ____ = [SELF.GLOBAL_num_dofs, *[ogm.GLOBAL_num_dofs for ogm in args]][:-1]
        upon = list()
        for j, _ in enumerate(args):
            upon.append(sum(____[0:j+1]))

        gvd = dict()
        for i in SELF:
            gvd_i_gv = SELF[i].full_vector
            for j, ogm in enumerate(args):
                gvd_i_gv = np.concatenate([gvd_i_gv, ogm[i].full_vector + upon[j]])
            gvd[i] = Gathering_Vector(i, gvd_i_gv)

        return SELF.__class__(gvd)

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = Gathering_Matrix_FIND(self._GM_)
        return self._find_