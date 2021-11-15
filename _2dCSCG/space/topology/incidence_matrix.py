import numpy as np
from SCREWS.frozen import FrozenOnly
from scipy.sparse import csc_matrix


class IncidenceMatrix(FrozenOnly):
    """
    Clearly, the incidence matrix of a form only depends on the local numbering
    and the basis function degree. We have fixed the way of numbering local
    dofs. That is the reason why we can evaluate basis here. Therefore, we can
    already compute the incidence matrix.

    """

    def __init__(self, FS):
        """ """
        assert FS.ndim == 2, " <IncidenceMatrix> "
        self._FS_ = FS
        self._freeze_self_()

    @property
    def _0Form_Inner(self):
        """grad.

        Here we generate the incidence matrix for inner 0-form in 2D.

        """
        sn = self._FS_.local_numbering._0Form_Inner
        dn = self._FS_.local_numbering._1Form_Inner
        E = np.zeros((self._FS_.num_basis._1Form_Inner[0],
                      self._FS_.num_basis._0Form_Inner[0]), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                E[dn[0][i,j], sn[0][i,j]]   = -1   # U
                E[dn[0][i,j], sn[0][i+1,j]] = +1   # D
        I, J = np.shape(dn[1])
        for j in range(J):
            for i in range(I):
                E[dn[1][i,j], sn[0][i,j]]   = -1    # L
                E[dn[1][i,j], sn[0][i,j+1]] = +1    # R
        return csc_matrix(E)

    @property
    def _1Form_Inner(self):
        """rot."""
        sn = self._FS_.local_numbering._1Form_Inner
        dn = self._FS_.local_numbering._2Form_Inner
        E = np.zeros((self._FS_.num_basis._2Form_Inner[0],
                      self._FS_.num_basis._1Form_Inner[0]), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                E[dn[0][i,j], sn[0][i,j  ]] = +1    # L
                E[dn[0][i,j], sn[0][i,j+1]] = -1    # R
                E[dn[0][i,j], sn[1][i  ,j]] = -1    # U
                E[dn[0][i,j], sn[1][i+1,j]] = +1    # D
        return csc_matrix(E)

    @property
    def _0Form_Outer(self):
        """curl.

        Here we generate the incidence matrix for inner 0-form in 2D.

        """
        sn = self._FS_.local_numbering._0Form_Outer
        dn = self._FS_.local_numbering._1Form_Outer
        E = np.zeros((self._FS_.num_basis._1Form_Outer[0],
                      self._FS_.num_basis._0Form_Outer[0]), dtype=int)
        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                E[dn[0][i,j], sn[0][i,j]]   = -1   # L ##########
                E[dn[0][i,j], sn[0][i,j+1]] = +1   # R ##########
        I, J = np.shape(dn[1])
        for j in range(J):
            for i in range(I):
                E[dn[1][i,j], sn[0][i,j]]   = +1    # U
                E[dn[1][i,j], sn[0][i+1,j]] = -1    # D
        return csc_matrix(E)

    @property
    def _1Form_Outer(self):
        """div."""
        sn = self._FS_.local_numbering._1Form_Outer
        dn = self._FS_.local_numbering._2Form_Outer
        E = np.zeros((self._FS_.num_basis._2Form_Outer[0],
                      self._FS_.num_basis._1Form_Outer[0]), dtype=int)

        I, J = np.shape(dn[0])
        for j in range(J):
            for i in range(I):
                E[dn[0][i,j], sn[1][i,j  ]] = -1    # L
                E[dn[0][i,j], sn[1][i,j+1]] = +1    # R
                E[dn[0][i,j], sn[0][i  ,j]] = -1    # U
                E[dn[0][i,j], sn[0][i+1,j]] = +1    # D
        return csc_matrix(E)