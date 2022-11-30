# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:12:47 2019
Aerodynamics, AE
TU Delft
"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from components.freeze.main import FrozenOnly
from scipy.sparse import csr_matrix


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
    def _3dCSCG_0Edge(self):
        """The selective matrix for the 0-edge-form."""
        p = self._FS_.p
        num_basis_0form = self._FS_.num_basis._3dCSCG_0Form[0]
        LNf = self._FS_.local_numbering._3dCSCG_0Form[0]
        LNe = self._FS_.local_numbering._3dCSCG_0Edge

        WB = np.zeros((p[0] + 1, num_basis_0form), dtype=int)
        WB[LNe['WB'], LNf[:,0,0].ravel('F')] = 1

        EB = np.zeros((p[0] + 1, num_basis_0form), dtype=int)
        EB[LNe['EB'], LNf[:,-1,0].ravel('F')] = 1

        WF = np.zeros((p[0] + 1, num_basis_0form), dtype=int)
        WF[LNe['WF'], LNf[:,0,-1].ravel('F')] = 1

        EF = np.zeros((p[0] + 1, num_basis_0form), dtype=int)
        EF[LNe['EF'], LNf[:,-1,-1].ravel('F')] = 1


        NB = np.zeros((p[1] + 1, num_basis_0form), dtype=int)
        NB[LNe['NB'], LNf[0, :, 0].ravel('F')] = 1

        SB = np.zeros((p[1] + 1, num_basis_0form), dtype=int)
        SB[LNe['SB'], LNf[-1, :, 0].ravel('F')] = 1

        NF = np.zeros((p[1] + 1, num_basis_0form), dtype=int)
        NF[LNe['NF'], LNf[0, :, -1].ravel('F')] = 1

        SF = np.zeros((p[1] + 1, num_basis_0form), dtype=int)
        SF[LNe['SF'], LNf[-1, :, -1].ravel('F')] = 1


        NW = np.zeros((p[2] + 1, num_basis_0form), dtype=int)
        NW[LNe['NW'], LNf[0, 0, :].ravel('F')] = 1

        SW = np.zeros((p[2] + 1, num_basis_0form), dtype=int)
        SW[LNe['SW'], LNf[-1, 0, :].ravel('F')] = 1

        NE = np.zeros((p[2] + 1, num_basis_0form), dtype=int)
        NE[LNe['NE'], LNf[0, -1, :].ravel('F')] = 1

        SE = np.zeros((p[2] + 1, num_basis_0form), dtype=int)
        SE[LNe['SE'], LNf[-1, -1, :].ravel('F')] = 1

        edge_wise_selective = {'WB': WB, 'EB': EB, 'WF': WF, 'EF': EF,
                               'NB': NB, 'SB': SB, 'NF': NF, 'SF': SF,
                               'NW': NW, 'SW': SW, 'NE': NE, 'SE': SE}

        mesh_element_wise_selective = np.vstack([WB, EB, WF, EF,
                                                 NB, SB, NF, SF,
                                                 NW, SW, NE, SE])
        mesh_element_wise_selective = csr_matrix(mesh_element_wise_selective)

        return mesh_element_wise_selective, edge_wise_selective

    @property
    def _3dCSCG_1Edge(self):
        """The selective matrix for the 1-edge-form."""
        p = self._FS_.p
        num_basis_1form = self._FS_.num_basis._3dCSCG_1Form[0]

        LNf = self._FS_.local_numbering._3dCSCG_1Form
        LNe = self._FS_.local_numbering._3dCSCG_1Edge


        WB = np.zeros((p[0], num_basis_1form), dtype=int)
        WB[LNe['WB'], LNf[0][:,0,0].ravel('F')] = 1

        EB = np.zeros((p[0], num_basis_1form), dtype=int)
        EB[LNe['EB'], LNf[0][:,-1,0].ravel('F')] = 1

        WF = np.zeros((p[0], num_basis_1form), dtype=int)
        WF[LNe['WF'], LNf[0][:,0,-1].ravel('F')] = 1

        EF = np.zeros((p[0], num_basis_1form), dtype=int)
        EF[LNe['EF'], LNf[0][:,-1,-1].ravel('F')] = 1


        NB = np.zeros((p[1], num_basis_1form), dtype=int)
        NB[LNe['NB'], LNf[1][0, :, 0].ravel('F')] = 1

        SB = np.zeros((p[1], num_basis_1form), dtype=int)
        SB[LNe['SB'], LNf[1][-1, :, 0].ravel('F')] = 1

        NF = np.zeros((p[1], num_basis_1form), dtype=int)
        NF[LNe['NF'], LNf[1][0, :, -1].ravel('F')] = 1

        SF = np.zeros((p[1], num_basis_1form), dtype=int)
        SF[LNe['SF'], LNf[1][-1, :, -1].ravel('F')] = 1


        NW = np.zeros((p[2], num_basis_1form), dtype=int)
        NW[LNe['NW'], LNf[2][0, 0, :].ravel('F')] = 1

        SW = np.zeros((p[2], num_basis_1form), dtype=int)
        SW[LNe['SW'], LNf[2][-1, 0, :].ravel('F')] = 1

        NE = np.zeros((p[2], num_basis_1form), dtype=int)
        NE[LNe['NE'], LNf[2][0, -1, :].ravel('F')] = 1

        SE = np.zeros((p[2], num_basis_1form), dtype=int)
        SE[LNe['SE'], LNf[2][-1, -1, :].ravel('F')] = 1

        edge_wise_selective = {'WB': WB, 'EB': EB, 'WF': WF, 'EF': EF,
                               'NB': NB, 'SB': SB, 'NF': NF, 'SF': SF,
                               'NW': NW, 'SW': SW, 'NE': NE, 'SE': SE}

        mesh_element_wise_selective = np.vstack([WB, EB, WF, EF,
                                                 NB, SB, NF, SF,
                                                 NW, SW, NE, SE])
        mesh_element_wise_selective = csr_matrix(mesh_element_wise_selective)

        return mesh_element_wise_selective, edge_wise_selective


    @property
    def _3dCSCG_0Trace(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._3dCSCG_0Form[0],
                       self._FS_.num_basis._3dCSCG_0Trace[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(
            self._FS_.num_basis._3dCSCG_0Trace[2])
        oln = self._FS_.local_numbering._3dCSCG_0Form[0]
        PU[oln[0, :, :].ravel('F'), sln['N']] = +1  # North
        PU[oln[-1, :, :].ravel('F'), sln['S']] = +1  # South
        PU[oln[:, 0, :].ravel('F'), sln['W']] = +1  # West
        PU[oln[:, -1, :].ravel('F'), sln['E']] = +1  # East
        PU[oln[:, :, 0].ravel('F'), sln['B']] = +1  # Back
        PU[oln[:, :, -1].ravel('F'), sln['F']] = +1  # Front
        N = PU.T
        nbt = self._FS_.num_basis._3dCSCG_0Trace[2]
        N_ = [None, None, None, None, None, None]
        m = 0
        for i in range(6):
            N_[i] = N[m:m + nbt['NSWEBF'[i]], :]
            m += nbt['NSWEBF'[i]]
        Nn, Ns, Nw, Ne, Nb, Nf = N_
        return csr_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne,
                               'B': Nb, 'F': Nf}

    @property
    def _3dCSCG_1Trace(self):
        """ """
        lnt = self._FS_.local_numbering._3dCSCG_1Trace
        lnf_dx, lnf_dy, lnf_dz = self._FS_.local_numbering._3dCSCG_1Form
        nbt = self._FS_.num_basis._3dCSCG_1Trace[2]
        nbf = self._FS_.num_basis._3dCSCG_1Form[0]
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
        return csr_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne, 'B': Nb, 'F': Nf}

    @property
    def _3dCSCG_2Trace(self):
        """ """
        PU = np.zeros((self._FS_.num_basis._3dCSCG_2Form[0],
                       self._FS_.num_basis._3dCSCG_2Trace[0]), dtype=int)
        sln = self.___generate_gathering_hybrid_element___(
            self._FS_.num_basis._3dCSCG_2Trace[2])
        oln_NS, oln_WE, oln_BF = self._FS_.local_numbering._3dCSCG_2Form
        PU[oln_NS[0, :, :].ravel('F'), sln['N']] = +1  # North
        PU[oln_NS[-1, :, :].ravel('F'), sln['S']] = +1  # South
        PU[oln_WE[:, 0, :].ravel('F'), sln['W']] = +1  # West
        PU[oln_WE[:, -1, :].ravel('F'), sln['E']] = +1  # East
        PU[oln_BF[:, :, 0].ravel('F'), sln['B']] = +1  # Back
        PU[oln_BF[:, :, -1].ravel('F'), sln['F']] = +1  # Front
        N = PU.T
        nbt = self._FS_.num_basis._3dCSCG_2Trace[2]
        N_ = [None, None, None, None, None, None]
        m = 0
        for i in range(6):
            N_[i] = N[m:m + nbt['NSWEBF'[i]], :]
            m += nbt['NSWEBF'[i]]
        Nn, Ns, Nw, Ne, Nb, Nf = N_
        return csr_matrix(N), {'N': Nn, 'S': Ns, 'W': Nw, 'E': Ne,
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




    @property
    def _3dCSCG_0LocalTrace(self):
        return self._3dCSCG_0Trace
    @property
    def _3dCSCG_1LocalTrace(self):
        return self._3dCSCG_1Trace
    @property
    def _3dCSCG_2LocalTrace(self):
        return self._3dCSCG_2Trace





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\space\topology\selective_matrix.py
    from objects.CSCG._3d.master import SpaceInvoker

    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])

    print(space._selective_matrix_._3dCSCG_1Edge)