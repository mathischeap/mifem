# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly



class SpaMat_Adjust(FrozenOnly):
    """Unlike `Customize`, `Adjust` will make the changes in real time and make a new
    EWC_SparseMatrix afterwards.
    """
    def __init__(self, spa_mat):
        """"""
        self._spa_mat_ = spa_mat
        self._freeze_self_()


    def clear_rows_according_to(
        self, itp, AS='local'):
        """We will clear rows of the sparse matrix M locally, and set
        M[pds.dofs, :] = 0

        Parameters
        ----------
        itp
        AS

        Returns
        -------

        """

        assert self._spa_mat_.elements._mesh_ == itp._mesh_, \
            "row itp elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if AS == 'local':

            LDF_row = itp.local.dofs

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_

                SMe = self._spa_mat_[e].copy().tolil()

                ROW = LDF_row[e]

                SMe[ROW, :] = 0

                SPA_MAT[e] = SMe.tocsr()

        else:
            raise NotImplementedError()

        for e in self._spa_mat_:
            if e not in SPA_MAT:
                SPA_MAT[e] = self._spa_mat_[e]

        SPA_MAT = self._spa_mat_.__class__(
            self._spa_mat_._elements_, SPA_MAT, cache_key_generator = 'no_cache')

        if self._spa_mat_._gathering_matrices_0_ is not None and \
            self._spa_mat_._gathering_matrices_1_ is not None:

            SPA_MAT.gathering_matrices = self._spa_mat_.gathering_matrices

        return SPA_MAT

    def identify_rows_according_to(
        self, itp, col_itp=None, AS='local'):
        """We will identify rows of the sparse matrix M locally, and set
        M[pds.dofs, :] = 0; M[pds.dofs, pds.dofs] = 1

        Note that it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        Parameters
        ----------
        itp
        col_itp
        AS

        Returns
        -------

        """
        if col_itp is not None:
            return self.___Pr_identify_rows_according_to___(itp, col_itp=col_itp, AS=AS)

        assert self._spa_mat_.elements._mesh_ == itp._mesh_, \
            "row PartialDofs elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if AS == 'local':

            LDF_row = itp.local.dofs

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_

                SMe = self._spa_mat_[e].copy().tolil()

                ROW = LDF_row[e]

                SMe[ROW, :] = 0
                SMe[ROW, ROW] = 1

                SPA_MAT[e] = SMe.tocsr()

        else:
            raise NotImplementedError()

        for e in self._spa_mat_:
            if e not in SPA_MAT:
                SPA_MAT[e] = self._spa_mat_[e]

        SPA_MAT = self._spa_mat_.__class__(
            self._spa_mat_._elements_, SPA_MAT, cache_key_generator = 'no_cache')

        if self._spa_mat_._gathering_matrices_0_ is not None and \
            self._spa_mat_._gathering_matrices_1_ is not None:
            SPA_MAT.gathering_matrices = self._spa_mat_.gathering_matrices

        return SPA_MAT

    def ___Pr_identify_rows_according_to___(
        self, row_itp, col_itp, AS='local'):
        """We will identify rows of the sparse matrix M locally, and set
        M[row_pds.dofs, :] = 0; M[row_pds.dofs, col_pds.dofs] = 1

        Note that it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        Parameters
        ----------
        row_itp
        col_itp
        AS

        Returns
        -------

        """

        assert self._spa_mat_.elements._mesh_ == row_itp._mesh_, \
            "row PartialDofs elements != EWC elements."
        assert self._spa_mat_.elements._mesh_ == col_itp._mesh_, \
            "col PartialDofs elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if AS == 'local':

            LDF_row = row_itp.local.dofs
            LDF_col = col_itp.local.dofs

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_ and e in LDF_col

                SMe = self._spa_mat_[e].copy().tolil()

                ROW = LDF_row[e]
                COL = LDF_col[e]

                SMe[ROW, :] = 0
                SMe[ROW, COL] = 1

                SPA_MAT[e] = SMe.tocsr()

        else:
            raise NotImplementedError()

        for e in self._spa_mat_:
            if e not in SPA_MAT:
                SPA_MAT[e] = self._spa_mat_[e]

        SPA_MAT = self._spa_mat_.__class__(
            self._spa_mat_._elements_, SPA_MAT, cache_key_generator = 'no_cache')

        if self._spa_mat_._gathering_matrices_0_ is not None and \
            self._spa_mat_._gathering_matrices_1_ is not None:
            SPA_MAT.gathering_matrices = self._spa_mat_.gathering_matrices

        return SPA_MAT

