# -*- coding: utf-8 -*-
import numpy as np
from components.freeze.main import FrozenOnly
from scipy.sparse import csc_matrix


class TraceMatrix(FrozenOnly):
    """This is the trace matrix.

    The N matrix will select from the cochain.local and get the cochain on the element boundary. And also the positive
    out-ward norm is included. So on the North, West, Back sides, we have '-', will on the 'South', 'East',
    'Front' sides, we will have '+'.

    IMPORTANT:
        Note here, the 'N'-matrix is for `mesh-element` not `trace-element`. So we already put 'N-north', 'N-south',
        'N-west', 'N-east', ...... together into a 'N'-matrix. And if the trace-form is not hybrid, then on the
        trace-element boundary, we should only have one dof, but that will be the case after assembling. Here, the
        assembling is not included. So for hybrid or non-hybrid trace-forms, we return the same N.

        Or in other words, for non-hybrid trace-form, the trace-boundary-dof have two (for 1-trace-form) or four
        (for 0-trace-form) rows referring to it.
    """

    def __init__(self, FS):
        """ """
        assert FS.ndim == 2, " <TraceMatrix> "
        self._FS_ = FS
        self._freeze_self_()

    @property
    def _2dCSCG_0Trace_Inner(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._2dCSCG_0Form_Inner[0], self._FS_.num_basis._0Trace_Inner[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(self._FS_.num_basis._0Trace_Inner[2])
        oln = self._FS_.local_numbering._2dCSCG_0Form_Inner[0]
        PU[oln[0, :].ravel('F'), sln['U']] = -1  # Upper
        PU[oln[-1, :].ravel('F'), sln['D']] = +1  # Down
        PU[oln[:, 0].ravel('F'), sln['L']] = -1  # Left
        PU[oln[:, -1].ravel('F'), sln['R']] = +1  # Right
        T = PU.T
        nbt = self._FS_.num_basis._0Trace_Inner[2]
        T_ = [None, None, None, None]
        m = 0
        for i in range(4):
            T_[i] = T[m:m + nbt['UDLR'[i]], :]
            m += nbt['UDLR'[i]]
        Tn, Ts, Tw, Te = T_
        return csc_matrix(T), {'U': Tn, 'D': Ts, 'L': Tw, 'R': Te}

    @property
    def _2dCSCG_0Trace_Outer(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._2dCSCG_0Form_Outer[0], self._FS_.num_basis._0Trace_Onner[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(self._FS_.num_basis._0Trace_Onner[2])
        oln = self._FS_.local_numbering._2dCSCG_0Form_Outer[0]
        PU[oln[0, :].ravel('F'), sln['U']] = +1  # Upper
        PU[oln[-1, :].ravel('F'), sln['D']] = -1  # Down
        PU[oln[:, 0].ravel('F'), sln['L']] = +1  # Left    ##########
        PU[oln[:, -1].ravel('F'), sln['R']] = -1  # Right   #############
        T = PU.T
        nbt = self._FS_.num_basis._0Trace_Onner[2]
        T_ = [None, None, None, None]
        m = 0
        for i in range(4):
            T_[i] = T[m:m + nbt['UDLR'[i]], :]
            m += nbt['UDLR'[i]]
        Tn, Ts, Tw, Te = T_
        return csc_matrix(T), {'U': Tn, 'D': Ts, 'L': Tw, 'R': Te}



    @property
    def _2dCSCG_1Trace_Inner(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._2dCSCG_1Form_Inner[0], self._FS_.num_basis._1Trace_Inner[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(self._FS_.num_basis._1Trace_Inner[2])
        oln_WE, oln_NS = self._FS_.local_numbering._2dCSCG_1Form_Inner
        PU[oln_NS[0, :].ravel('F'), sln['U']] = -1  # Upper
        PU[oln_NS[-1, :].ravel('F'), sln['D']] = +1  # Down
        PU[oln_WE[:, 0].ravel('F'), sln['L']] = +1  # Left
        PU[oln_WE[:, -1].ravel('F'), sln['R']] = -1  # Right
        T = PU.T
        nbt = self._FS_.num_basis._1Trace_Inner[2]
        T_ = [None, None, None, None]
        m = 0
        for i in range(4):
            T_[i] = T[m:m + nbt['UDLR'[i]], :]
            m += nbt['UDLR'[i]]
        Tn, Ts, Tw, Te = T_
        return csc_matrix(T), {'U': Tn, 'D': Ts, 'L': Tw, 'R': Te}

    @property
    def _2dCSCG_1Trace_Outer(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._2dCSCG_1Form_Outer[0], self._FS_.num_basis._2dCSCG_1Trace_Outer[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(self._FS_.num_basis._2dCSCG_1Trace_Outer[2])
        oln_NS, oln_WE = self._FS_.local_numbering._2dCSCG_1Form_Outer
        PU[oln_NS[ 0, :].ravel('F'), sln['U']] = +1  # Upper
        PU[oln_NS[-1, :].ravel('F'), sln['D']] = -1  # Down
        PU[oln_WE[:,  0].ravel('F'), sln['L']] = -1  # Left
        PU[oln_WE[:, -1].ravel('F'), sln['R']] = +1  # Right


        T = PU.T
        nbt = self._FS_.num_basis._2dCSCG_1Trace_Outer[2]
        T_ = [None, None, None, None]
        m = 0
        for i in range(4):
            T_[i] = T[m:m + nbt['UDLR'[i]], :]
            m += nbt['UDLR'[i]]

        TU, TD, TL, TR = T_
        return csc_matrix(T), {'U': TU, 'D': TD, 'L': TL, 'R': TR}



    @staticmethod
    def ___generate_gathering_hybrid_element___(basis_onsides):
        """ """
        si2n = ('U', 'D', 'L', 'R')  # the edge numbering sequence: U->D->L->R
        _gte_ = {}
        current_num = 0
        for sn in si2n:
            _gte_[sn] = np.arange(current_num, current_num + basis_onsides[sn]).tolist()
            current_num += basis_onsides[sn]
        return _gte_
