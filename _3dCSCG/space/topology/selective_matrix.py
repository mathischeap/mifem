# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:12:47 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from SCREWS.frozen import FrozenOnly
from scipy.sparse import csc_matrix


class SelectiveMatrix(FrozenOnly):
    """
    This is the selective matrix. It is like the Trace_matrix, but without '-'.

    """

    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <SelectiveMatrix> "
        self._FS_ = FS
        self._freeze_self_()

    @property
    def _0Trace(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._0Form[0],
                       self._FS_.num_basis._0Trace[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(
            self._FS_.num_basis._0Trace[2])
        oln = self._FS_.local_numbering._0Form[0]
        PU[oln[0, :, :].ravel('F'), sln['N']] = +1  # North
        PU[oln[-1, :, :].ravel('F'), sln['S']] = +1  # South
        PU[oln[:, 0, :].ravel('F'), sln['W']] = +1  # West
        PU[oln[:, -1, :].ravel('F'), sln['E']] = +1  # East
        PU[oln[:, :, 0].ravel('F'), sln['B']] = +1  # Back
        PU[oln[:, :, -1].ravel('F'), sln['F']] = +1  # Front
        N = PU.T
        nbt = self._FS_.num_basis._0Trace[2]
        N_ = [None, None, None, None, None, None]
        m = 0
        for i in range(6):
            N_[i] = N[m:m + nbt['NSWEBF'[i]], :]
            m += nbt['NSWEBF'[i]]
        Nn, Ns, Nw, Ne, Nb, Nf = N_
        return csc_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne,
                               'B': Nb, 'F': Nf}

    @property
    def _1Trace(self):
        """ """
        lnt = self._FS_.local_numbering._1Trace
        lnf_dx, lnf_dy, lnf_dz = self._FS_.local_numbering._1Form
        nbt = self._FS_.num_basis._1Trace[2]
        nbf = self._FS_.num_basis._1Form[0]
        # To see why we have following +,-, see the figure: IMG_20190501_170517.jpg in
        # folder .\bin\miscellaneous
        # ____ North __________________________________________________________________
        Nn = np.zeros((nbt['N'], nbf), dtype=int)
        Nn[lnt['N'][0].ravel('F'), lnf_dy[0, :, :].ravel('F')] = +1
        Nn[lnt['N'][1].ravel('F'), lnf_dz[0, :, :].ravel('F')] = +1
        # ____ South __________________________________________________________________
        Ns = np.zeros((nbt['S'], nbf), dtype=int)
        Ns[lnt['S'][0].ravel('F'), lnf_dy[-1, :, :].ravel('F')] = +1
        Ns[lnt['S'][1].ravel('F'), lnf_dz[-1, :, :].ravel('F')] = +1
        # ____ West ___________________________________________________________________
        Nw = np.zeros((nbt['W'], nbf), dtype=int)
        Nw[lnt['W'][0].ravel('F'), lnf_dx[:, 0, :].ravel('F')] = +1
        Nw[lnt['W'][1].ravel('F'), lnf_dz[:, 0, :].ravel('F')] = +1
        # ____ East ___________________________________________________________________
        Ne = np.zeros((nbt['E'], nbf), dtype=int)
        Ne[lnt['E'][0].ravel('F'), lnf_dx[:, -1, :].ravel('F')] = +1
        Ne[lnt['E'][1].ravel('F'), lnf_dz[:, -1, :].ravel('F')] = +1
        # ____ Back ___________________________________________________________________
        Nb = np.zeros((nbt['B'], nbf), dtype=int)
        Nb[lnt['B'][0].ravel('F'), lnf_dx[:, :, 0].ravel('F')] = +1
        Nb[lnt['B'][1].ravel('F'), lnf_dy[:, :, 0].ravel('F')] = +1
        # ____ Front __________________________________________________________________
        Nf = np.zeros((nbt['F'], nbf), dtype=int)
        Nf[lnt['F'][0].ravel('F'), lnf_dx[:, :, -1].ravel('F')] = +1
        Nf[lnt['F'][1].ravel('F'), lnf_dy[:, :, -1].ravel('F')] = +1
        # ------------------------------------------------------------------------------
        N = np.vstack((Nn, Ns, Nw, Ne, Nb, Nf))
        return csc_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne, 'B': Nb, 'F': Nf}

    @property
    def _2Trace(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._2Form[0],
                       self._FS_.num_basis._2Trace[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(
            self._FS_.num_basis._2Trace[2])
        oln_NS, oln_WE, oln_BF = self._FS_.local_numbering._2Form
        PU[oln_NS[0, :, :].ravel('F'), sln['N']] = +1  # North
        PU[oln_NS[-1, :, :].ravel('F'), sln['S']] = +1  # South
        PU[oln_WE[:, 0, :].ravel('F'), sln['W']] = +1  # West
        PU[oln_WE[:, -1, :].ravel('F'), sln['E']] = +1  # East
        PU[oln_BF[:, :, 0].ravel('F'), sln['B']] = +1  # Back
        PU[oln_BF[:, :, -1].ravel('F'), sln['F']] = +1  # Front
        N = PU.T
        nbt = self._FS_.num_basis._2Trace[2]
        N_ = [None, None, None, None, None, None]
        m = 0
        for i in range(6):
            N_[i] = N[m:m + nbt['NSWEBF'[i]], :]
            m += nbt['NSWEBF'[i]]
        Nn, Ns, Nw, Ne, Nb, Nf = N_
        return csc_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne,
                               'B': Nb, 'F': Nf}

    @staticmethod
    def ___generate_gathering_hybrid_element___(basis_on_sides):
        """ """
        si2n = ('N', 'S', 'W', 'E', 'B', 'F')
        _gte_ = {}
        current_num = 0
        for sn in si2n:
            _gte_[sn] = np.arange(current_num,
                                  current_num + basis_on_sides[
                                      sn]).tolist()
            current_num += basis_on_sides[sn]
        return _gte_