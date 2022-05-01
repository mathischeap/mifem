


from screws.decorators.accepts import accepts
from screws.freeze.base import FrozenOnly



class SpaMat_Adjust(FrozenOnly):
    """Unlike `Customize`, `Adjust` will make the changes in real time and make a new
    EWC_SparseMatrix immediately.
    """
    def __init__(self, spa_mat):
        """"""
        self._spa_mat_ = spa_mat
        self._freeze_self_()


    @accepts('self',
             ('PartialDofs', 'PartialCochain'),
             ('PartialDofs', 'PartialCochain'))
    def identify_rows_according_to_two_CSCG_partial_dofs(
        self, row_pds, col_pds, interpreted_as='local_dofs'):
        """We will identify rows of the sparse matrix M locally, and set
        M[row_pds.dofs, :] = 0; M[row_pds.dofs, col_pds.dofs] = 1

        Note that it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        Parameters
        ----------
        row_pds
        col_pds
        interpreted_as

        Returns
        -------

        """
        if row_pds.__class__.__name__ == 'PartialCochain':
            row_pds = row_pds.dofs
        if col_pds.__class__.__name__ == 'PartialCochain':
            col_pds = col_pds.dofs

        assert self._spa_mat_.elements._mesh_ == row_pds._mesh_, \
            "row PartialDofs elements != EWC elements."
        assert self._spa_mat_.elements._mesh_ == col_pds._mesh_, \
            "col PartialDofs elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if interpreted_as == 'local_dofs':

            LDF_row = row_pds.interpreted_as.local_dofs
            LDF_col = col_pds.interpreted_as.local_dofs

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


    @accepts('self', ('PartialDofs', 'PartialCochain'))
    def clear_rows_according_to_CSCG_partial_dofs(
        self, pds, interpreted_as='local_dofs'):
        """We will clear rows of the sparse matrix M locally, and set
        M[pds.dofs, :] = 0

        Parameters
        ----------
        pds
        interpreted_as

        Returns
        -------

        """
        if pds.__class__.__name__ == 'PartialCochain':
            pds = pds.dofs

        assert self._spa_mat_.elements._mesh_ == pds._mesh_, \
            "row PartialDofs elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if interpreted_as == 'local_dofs':

            LDF_row = pds.interpreted_as.local_dofs

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

    @accepts('self', ('PartialDofs', 'PartialCochain'))
    def identify_rows_according_to_CSCG_partial_dofs(
        self, pds, interpreted_as='local_dofs'):
        """We will identify rows of the sparse matrix M locally, and set
        M[pds.dofs, :] = 0; M[pds.dofs, pds.dofs] = 1

        Note that it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        Parameters
        ----------
        pds
        interpreted_as

        Returns
        -------

        """
        if pds.__class__.__name__ == 'PartialCochain':
            pds = pds.dofs

        assert self._spa_mat_.elements._mesh_ == pds._mesh_, \
            "row PartialDofs elements != EWC elements."

        SPA_MAT = dict()  # the data of the new EWC_SparseMatrix.

        if interpreted_as == 'local_dofs':

            LDF_row = pds.interpreted_as.local_dofs

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