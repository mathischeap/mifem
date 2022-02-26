# -*- coding: utf-8 -*-
"""
"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from screws.frozen import FrozenOnly
from importlib import import_module
from tools.deprecated.linear_system.lhs_rhs import LHS, RHS



class LinearSystem(FrozenOnly):
    """
    A linear system.

    :param lhs:
    :param rhs:
    """
    def __init__(self, lhs, rhs=None, name='unnamed linear system'):
        lhs = list(lhs) if isinstance(lhs, tuple) else lhs
        assert isinstance(lhs, list), "lhs has to be list"

        LEN = list()
        for i, lhs_i in enumerate(lhs):
            LEN.append(len(lhs_i))
            lhs[i] = list(lhs_i)

        for L in LEN:
            assert L == LEN[0], "I need a 2d lhs!"


        lhs = tuple(lhs)
        s0 = len(lhs)
        s1 = len(lhs[0])
        self._shape_ = (s0, s1)

        # self._shape_ = np.shape(lhs)
        # print(self._shape_ )

        if rhs is None:
            # in this case, we actually do not want to solve the system. We will just want to have
            # the lhs and study its like eigen values, eigen vectors and condition number.
            rhs = list()
            for i in range(self._shape_[0]):
                rhs.append(None)
        else:
            rhs = list(rhs) if isinstance(rhs, tuple) else rhs
            assert isinstance(rhs, list), " rhs must be a list."
            assert len(rhs) == self._shape_[0], "rhs shape does not match lhs shape."

        if self._shape_ == (1,1):
            assert lhs[0][0] is not None

        self._lhs_ = LHS(self, lhs)
        self._rhs_ = RHS(self, rhs)
        self.name = name
        self._solve_ = _LinearSystem_Solve(self)
        self._structure_ = _LinearSystem_Structure(self)
        self._condition_ = _LinearSystem_Condition(self)
        self._visualize_ = _LinearSystem_Visualize(self)
        self._freeze_self_()

    @property
    def name(self):
        return self._name_

    @name.setter
    def name(self, name):
        assert isinstance(name, str), "name must be str."
        self._name_ = name

    @property
    def lhs(self):
        return self._lhs_

    @property
    def rhs(self):
        return self._rhs_

    @property
    def shape(self):
        return self._shape_



    @property
    def structure(self):
        return self._structure_

    @property
    def solve(self):
        return self._solve_

    @property
    def condition(self):
        return self._condition_

    @property
    def visualize(self):
        return self._visualize_



class _LinearSystem_Solve(FrozenOnly):
    def __init__(self, LS):
        self._LS_ = LS
        self._solver_instance_ = None
        self._results_ = None
        self._message_ = None
        self._info_ = None
        self._others_ = None
        self._freeze_self_()

    def __call__(self, solver_name, routine_name=None):
        """

        :param solver_name: The name of the solver we use.
        :param routine_name: The routine we choose.
        :return:
        """
        if self._solver_instance_ is None or self._solver_instance_.__class__.__name__ != solver_name:
            self._solver_instance_ = getattr(
                import_module(
                    'tools.deprecated.linear_system.solvers'), solver_name)(self._LS_)
        else:
            pass
        if routine_name is None:
            routine_name = self._solver_instance_._default_routine_
        return getattr(self._solver_instance_, routine_name)

    @property
    def results(self):
        return self._results_

    @property
    def message(self):
        return self._message_

    @property
    def info(self):
        return self._info_

    @property
    def others(self):
        return self._others_


class _LinearSystem_Condition(FrozenOnly):
    def __init__(self, LS):
        self._LS_ = LS
        self._cn_ = None
        self._freeze_self_()

    @property
    def condition_number(self):
        """"""
        assert self._LS_.structure.IS_all_global_matrices, "Need IS_all_global_matrices."
        if self._cn_ is None:
            A = self._LS_.lhs.___PRIVATE_convert_into_one_sparse_A___()
            M = A.___PRIVATE_gather_M_to_core___()
            if rAnk == mAster_rank:
                M = M.toarray()
                self._cn_ = np.linalg.cond(M)
            else:
                self._cn_ = None
            self._cn_ = cOmm.bcast(self._cn_, root=mAster_rank)
        return self._cn_


class _LinearSystem_Structure(FrozenOnly):
    def __init__(self, LS):
        self._LS_ = LS
        self._freeze_self_()

    @property
    def regularity(self):
        """
        (str) Return a key that represents the structure of the linear system.

        This must be real time, as we may change the regularity time by time.
        """
        lhs_structure = self._LS_.lhs.DO_parse_structure()
        rhs_structure = self._LS_.rhs.DO_parse_structure()
        regularity = ''
        for i in range(self._LS_.shape[0]):
            regularity += ' '.join(lhs_structure[i])
            regularity += ' | '
            regularity += rhs_structure[i]
            regularity += '\n'
        return regularity

    @property
    def regularity_ravel(self):
        regularity = self.regularity
        regularity = regularity.replace(' ', '')
        regularity = regularity.replace('\n', '')
        return regularity

    @property
    def IS_all_global_matrices(self):
        """(bool) Return ``True`` if all left hand side blocks are global matrices."""
        regularity1 = self._LS_.lhs.___PRIVATE_IS_all_global_matrices___
        regularity2 = self._LS_.rhs.___PRIVATE_IS_all_global_matrices___
        return regularity1 and regularity2



class _LinearSystem_Visualize(FrozenOnly):
    def __init__(self, LS):
        self._LS_ = LS
        self._freeze_self_()

    def spy(self, **kwargs):
        A = self._LS_.lhs.___PRIVATE_convert_into_one_sparse_A___()
        A.visualize.spy(**kwargs)





if __name__ == '__main__':
    # mpiexec python TOOLS\linear_system\main.py
    lhs00 = 0
    lhs01 = 1
    lhs10 = 2
    lhs11 = 3

    lhs = [[lhs00, lhs01],
           [lhs10, lhs11]]
    A = LinearSystem(lhs)
