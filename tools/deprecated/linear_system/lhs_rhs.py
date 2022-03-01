# -*- coding: utf-8 -*-
from screws.frozen import FrozenOnly
from scipy import sparse as spspa
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector, GlobalMatrix




class LHS(FrozenOnly):
    def __init__(self, LS, lhs):
        self._LS_ = LS
        self._lhs_ = lhs
        self._shape_ = LS._shape_
        self.DO_parse_structure()
        self._freeze_self_()

    def __getitem__(self, item):
        return self._lhs_[item]


    def DO_parse_structure(self):
        """
        This is very important.

        Whenever this linear system class support new data structure, we need
        to update this method.

        :return: The current structure.
        :rtype: list
        """
        structure = [['------------' for _ in range(self._shape_[1])] for _ in range(self._shape_[0])]

        s0 = len(self._lhs_)
        s1 = len(self._lhs_[0])
        assert self._shape_ == (s0, s1), "shape changed?"

        for i in range(self._shape_[0]):
            for j in range(self._shape_[1]):
                lhs_ij = self._lhs_[i][j]
                if lhs_ij.__class__.__name__ == 'GlobalMatrix':
                    structure[i][j] = '-----GM-----'
                elif lhs_ij.__class__.__name__ == 'EWC_SparseMatrix':
                    structure[i][j] = 'EWC_SpaseMat'
                elif lhs_ij is None:
                    pass
                else:
                    raise Exception(f'Linear system lhs[{i}][{j}] is {lhs_ij}, wrong.')

        return structure


    def ___PRIVATE_convert_into_one_sparse_A___(self):
        """"""
        assert self.___PRIVATE_IS_all_global_matrices___
        BMAT = [[1 for _ in range(self._shape_[1])] for _ in range(self._shape_[0])]
        for i in range(self._shape_[0]):
            for j in range(self._shape_[1]):
                if self._lhs_[i][j] is None:
                    # noinspection PyTypeChecker
                    BMAT[i][j] = None
                else:
                    BMAT[i][j] = self._lhs_[i][j].M
        return GlobalMatrix(spspa.bmat(BMAT, format='csc'))

    @property
    def ___PRIVATE_IS_all_global_matrices___(self):
        """(bool) Return ``True`` if all left hand side blocks are global matrices."""
        regularity = self.DO_parse_structure()
        for ri in regularity:
            for rij in ri:
                if rij in ('-----GM-----', '------------'):
                    pass
                else:
                    return False
        return True





#--------------------------------------------------------------------------------------------------------


class RHS(FrozenOnly):
    def __init__(self, LS, rhs):
        self._LS_ = LS
        self._rhs_ = rhs
        self.DO_parse_structure()
        self._freeze_self_()

    def __getitem__(self, item):
        return self._rhs_[item]

    def __setitem__(self, key, value):
        self._rhs_[key] = value

    def DO_parse_structure(self):
        """
        This is very important.

        Whenever this linear system class support new data structure, we need
        to update this method.

        :return: The current structure.
        :rtype: list
        """
        structure = list()
        assert len(self._rhs_) == self._LS_.shape[0]
        for i in range(self._LS_.shape[0]):
            if self._rhs_[i].__class__.__name__ == 'GlobalVector':
                structure.append('gv')
            elif self._rhs_[i].__class__.__name__ == 'EWC_ColumnVector':
                structure.append('EWC_ColVec')
            elif self._rhs_[i] is None:
                structure.append('-')
            else:
                raise Exception(f'Linear system rhs[{i}] is {self._rhs_[i]}, wrong.')
        return structure





    def ___PRIVATE_convert_into_one_sparse_b___(self):
        """"""
        assert self.___PRIVATE_IS_all_global_matrices___
        VSTACK = [self._rhs_[i].V for i in range(len(self._rhs_))]
        return GlobalVector(spspa.vstack(VSTACK, format='csc'))

    @property
    def ___PRIVATE_IS_all_global_matrices___(self):
        """(bool) Return ``True`` if all left hand side blocks are global matrices."""
        regularity = self.DO_parse_structure()
        for ri in regularity:
            if ri == 'gv': # wo do not take None as global cause we do not now the shape.
                pass
            else:
                return False
        return True