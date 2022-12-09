# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:50 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from components.freeze.main import FrozenOnly

class LocalNumbering(FrozenOnly):
    """ """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <LocalNumbering> "
        self._FS_ = FS
        self._cache_E0_ = dict()
        self._cache_E1_ = dict()
        self._cache_T0_ = dict()
        self._cache_T1_ = dict()
        self._freeze_self_()






    @property
    def _3dCSCG_0Edge(self):
        """Return the edge-element-wise (NOT mesh-element-wise) local numbering of 0-edge-form."""
        p = self._FS_.p
        NS = [_ for _ in range(p[0]+1)]
        WE = [_ for _ in range(p[1]+1)]
        BF = [_ for _ in range(p[2]+1)]
        EEW__LN = \
            {'WB': NS,
             'EB': NS,
             'WF': NS,
             'EF': NS,
             'NB': WE,
             'SB': WE,
             'NF': WE,
             'SF': WE,
             'NW': BF,
             'SW': BF,
             'NE': BF,
             'SE': BF}
        return EEW__LN

    @property
    def _3dCSCG_1Edge(self):
        """Return the edge-element-wise (NOT mesh-element-wise) local numbering of 1-edge-form."""
        p = self._FS_.p
        NS = [_ for _ in range(p[0])]
        WE = [_ for _ in range(p[1])]
        BF = [_ for _ in range(p[2])]
        EEW__LN = \
            {'WB': NS,
             'EB': NS,
             'WF': NS,
             'EF': NS,
             'NB': WE,
             'SB': WE,
             'NF': WE,
             'SF': WE,
             'NW': BF,
             'SW': BF,
             'NE': BF,
             'SE': BF}
        return EEW__LN

    def ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_0edge_edge___(self, edge_name):
        """"""
        if edge_name not in self._cache_E0_:

            px, py, pz = self._FS_.p
            px += 1
            py += 1
            pz += 1

            if edge_name == 'WB':
                DOFS = np.arange(px)
            elif edge_name == 'EB':
                DOFS = np.arange(px) + px
            elif edge_name == 'WF':
                DOFS = np.arange(px) + px * 2
            elif edge_name == 'EF':
                DOFS = np.arange(px) + px * 3
            elif edge_name == 'NB':
                DOFS = np.arange(py) + px * 4
            elif edge_name == 'SB':
                DOFS = np.arange(py) + px * 4 + py
            elif edge_name == 'NF':
                DOFS = np.arange(py) + px * 4 + py * 2
            elif edge_name == 'SF':
                DOFS = np.arange(py) + px * 4 + py * 3
            elif edge_name == 'NW':
                DOFS = np.arange(pz) + px * 4 + py * 4
            elif edge_name == 'SW':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz
            elif edge_name == 'NE':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz * 2
            elif edge_name == 'SE':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz * 3
            else:
                raise Exception()

            self._cache_E0_[edge_name] = DOFS

        return self._cache_E0_[edge_name]

    def ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_1edge_edge___(self, edge_name):
        """"""
        if edge_name not in self._cache_E1_:

            px, py, pz = self._FS_.p

            if edge_name == 'WB':
                DOFS = np.arange(px)
            elif edge_name == 'EB':
                DOFS = np.arange(px) + px
            elif edge_name == 'WF':
                DOFS = np.arange(px) + px * 2
            elif edge_name == 'EF':
                DOFS = np.arange(px) + px * 3
            elif edge_name == 'NB':
                DOFS = np.arange(py) + px * 4
            elif edge_name == 'SB':
                DOFS = np.arange(py) + px * 4 + py
            elif edge_name == 'NF':
                DOFS = np.arange(py) + px * 4 + py * 2
            elif edge_name == 'SF':
                DOFS = np.arange(py) + px * 4 + py * 3
            elif edge_name == 'NW':
                DOFS = np.arange(pz) + px * 4 + py * 4
            elif edge_name == 'SW':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz
            elif edge_name == 'NE':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz * 2
            elif edge_name == 'SE':
                DOFS = np.arange(pz) + px * 4 + py * 4 + pz * 3
            else:
                raise Exception()

            self._cache_E1_[edge_name] = DOFS

        return self._cache_E1_[edge_name]


    @property
    def _3dCSCG_0Trace(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 0-trace-form."""
        p = self._FS_.p
        TEW__LN = {'N': (np.arange((p[1]+1)*(p[2]+1)).reshape((p[1]+1, p[2]+1), order='F'),),
                   'S': (np.arange((p[1]+1)*(p[2]+1)).reshape((p[1]+1, p[2]+1), order='F'),),
                   'W': (np.arange((p[0]+1)*(p[2]+1)).reshape((p[0]+1, p[2]+1), order='F'),),
                   'E': (np.arange((p[0]+1)*(p[2]+1)).reshape((p[0]+1, p[2]+1), order='F'),),
                   'B': (np.arange((p[0]+1)*(p[1]+1)).reshape((p[0]+1, p[1]+1), order='F'),),
                   'F': (np.arange((p[0]+1)*(p[1]+1)).reshape((p[0]+1, p[1]+1), order='F'),)}
        return TEW__LN
    
    @property
    def _3dCSCG_1Trace(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 1-trace-form."""
        p = self._FS_.p
        TEW__LN = dict()
        TEW__LN['N'] = (np.arange(p[1]*(p[2]+1)).reshape((p[1], p[2]+1), order='F'),
                        np.arange((p[1]+1)*p[2]).reshape((p[1]+1, p[2]), order='F') + p[1]*(p[2]+1))
        TEW__LN['S'] = (np.arange(p[1]*(p[2]+1)).reshape((p[1], p[2]+1), order='F'),
                        np.arange((p[1]+1)*p[2]).reshape((p[1]+1, p[2]), order='F') + p[1]*(p[2]+1))
        TEW__LN['W'] = (np.arange(p[0]*(p[2]+1)).reshape((p[0], p[2]+1), order='F'),
                        np.arange((p[0]+1)*p[2]).reshape((p[0]+1, p[2]), order='F') + p[0]*(p[2]+1))
        TEW__LN['E'] = (np.arange(p[0]*(p[2]+1)).reshape((p[0], p[2]+1), order='F'),
                        np.arange((p[0]+1)*p[2]).reshape((p[0]+1, p[2]), order='F') + p[0]*(p[2]+1))
        TEW__LN['B'] = (np.arange(p[0]*(p[1]+1)).reshape((p[0], p[1]+1), order='F'),
                        np.arange((p[0]+1)*p[1]).reshape((p[0]+1, p[1]), order='F') + p[0]*(p[1]+1))
        TEW__LN['F'] = (np.arange(p[0]*(p[1]+1)).reshape((p[0], p[1]+1), order='F'),
                        np.arange((p[0]+1)*p[1]).reshape((p[0]+1, p[1]), order='F') + p[0]*(p[1]+1))
        return TEW__LN
    
    @property
    def _3dCSCG_2Trace(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 2-trace-form."""
        p = self._FS_.p
        TEW__LN = {'N': (np.arange(p[1]*p[2]).reshape((p[1], p[2]), order='F'),),
                   'S': (np.arange(p[1]*p[2]).reshape((p[1], p[2]), order='F'),),
                   'W': (np.arange(p[0]*p[2]).reshape((p[0], p[2]), order='F'),),
                   'E': (np.arange(p[0]*p[2]).reshape((p[0], p[2]), order='F'),),
                   'B': (np.arange(p[0]*p[1]).reshape((p[0], p[1]), order='F'),),
                   'F': (np.arange(p[0]*p[1]).reshape((p[0], p[1]), order='F'),)}
        return TEW__LN




    def ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_0Trace_edge___(
        self, trace_face, trace_edge, hybrid=True):
        """We try to find the local dofs of 0-trace-form on the `trace_edge` of a trace-element on
        the  `trace_face` side of mesh-element.

        Parameters
        ----------
        trace_face
        trace_edge
        hybrid:
            If the trace form is hybrid?

        Returns
        -------

        """
        if hybrid:
            TFE = trace_face + trace_edge
            if TFE not in self._cache_T0_:

                TEW_LN = self._3dCSCG_0Trace

                LN = TEW_LN[trace_face]

                if trace_face in ('N', 'S'):
                    if trace_edge == 'W':
                        TEW_dofs = LN[0][0,:]
                    elif trace_edge == 'E':
                        TEW_dofs = LN[0][-1,:]
                    elif trace_edge == 'B':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'F':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                elif trace_face in ('W', 'E'):
                    if trace_edge == 'N':
                        TEW_dofs = LN[0][0,:]
                    elif trace_edge == 'S':
                        TEW_dofs = LN[0][-1,:]
                    elif trace_edge == 'B':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'F':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                elif trace_face in ('B', 'F'):
                    if trace_edge == 'N':
                        TEW_dofs = LN[0][0,:]
                    elif trace_edge == 'S':
                        TEW_dofs = LN[0][-1,:]
                    elif trace_edge == 'W':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'E':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                else:
                    raise Exception()

                num_basis_onsides = self._FS_.num_basis._3dCSCG_0Trace[2]
                num_basis_NS = num_basis_onsides['N']
                num_basis_WE = num_basis_onsides['W']
                num_basis_BF = num_basis_onsides['B']

                if trace_face =='N':
                    self._cache_T0_[TFE] = TEW_dofs
                elif trace_face == 'S':
                    self._cache_T0_[TFE] = TEW_dofs + num_basis_NS
                elif trace_face == 'W':
                    self._cache_T0_[TFE] = TEW_dofs + num_basis_NS * 2
                elif trace_face == 'E':
                    self._cache_T0_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE
                elif trace_face == 'B':
                    self._cache_T0_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE * 2
                elif trace_face == 'F':
                    self._cache_T0_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE * 2 + num_basis_BF
                else:
                    raise Exception()

            return self._cache_T0_[TFE]

        else:
            raise NotImplementedError()

    def ___PRIVATE_find_MESH_ELEMENT_WISE_local_dofs_of_1Trace_edge___(
        self, trace_face, trace_edge, hybrid=True):
        """We try to find the local dofs of 1-trace-form on the `trace_edge` of a trace-element on
        the  `trace_face` side of mesh-element.

        Parameters
        ----------
        trace_face
        trace_edge
        hybrid:
            If the trace form is hybrid?

        Returns
        -------

        """
        if hybrid:
            TFE = trace_face + trace_edge
            if TFE not in self._cache_T1_:

                TEW_LN = self._3dCSCG_1Trace

                LN = TEW_LN[trace_face]

                if trace_face in ('N', 'S'):
                    if trace_edge == 'W':
                        TEW_dofs = LN[1][0,:]
                    elif trace_edge == 'E':
                        TEW_dofs = LN[1][-1,:]
                    elif trace_edge == 'B':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'F':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                elif trace_face in ('W', 'E'):
                    if trace_edge == 'N':
                        TEW_dofs = LN[1][0,:]
                    elif trace_edge == 'S':
                        TEW_dofs = LN[1][-1,:]
                    elif trace_edge == 'B':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'F':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                elif trace_face in ('B', 'F'):
                    if trace_edge == 'N':
                        TEW_dofs = LN[1][0,:]
                    elif trace_edge == 'S':
                        TEW_dofs = LN[1][-1,:]
                    elif trace_edge == 'W':
                        TEW_dofs = LN[0][:,0]
                    elif trace_edge == 'E':
                        TEW_dofs = LN[0][:,-1]
                    else:
                        raise Exception()

                else:
                    raise Exception()

                num_basis_onsides = self._FS_.num_basis._3dCSCG_1Trace[2]
                num_basis_NS = num_basis_onsides['N']
                num_basis_WE = num_basis_onsides['W']
                num_basis_BF = num_basis_onsides['B']

                if trace_face =='N':
                    self._cache_T1_[TFE] = TEW_dofs
                elif trace_face == 'S':
                    self._cache_T1_[TFE] = TEW_dofs + num_basis_NS
                elif trace_face == 'W':
                    self._cache_T1_[TFE] = TEW_dofs + num_basis_NS * 2
                elif trace_face == 'E':
                    self._cache_T1_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE
                elif trace_face == 'B':
                    self._cache_T1_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE * 2
                elif trace_face == 'F':
                    self._cache_T1_[TFE] = TEW_dofs + num_basis_NS * 2 + num_basis_WE * 2 + num_basis_BF
                else:
                    raise Exception()

            return self._cache_T1_[TFE]

        else:
            raise NotImplementedError()




    @property
    def _3dCSCG_0LocalTrace(self):
        return self._3dCSCG_0Trace
    @property
    def _3dCSCG_1LocalTrace(self):
        return self._3dCSCG_1Trace
    @property
    def _3dCSCG_2LocalTrace(self):
        return self._3dCSCG_2Trace









    @property
    def _3dCSCG_0Form(self):
        """Return the mesh-element-wise local numbering of 0-form."""
        p = [self._FS_.p[i]+1 for i in range(self._FS_.ndim)]
        _ln_ = (np.arange(self._FS_.num_basis._3dCSCG_0Form[0]).reshape(*p, order='F'),)
        return _ln_
    
    @property
    def _3dCSCG_1Form(self):
        """Return the mesh-element-wise local numbering of 1-form."""
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j]+1 for j in range(self._FS_.ndim)]
            p[i] -= 1
            I = 0 if i ==0 else np.sum(self._FS_.num_basis._3dCSCG_1Form[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._3dCSCG_1Form[1][i]).reshape(p, order='F') + I,)
        return _ln_
    
    @property
    def _3dCSCG_2Form(self):
        """Return the mesh-element-wise local numbering of 2-form."""
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j] for j in range(self._FS_.ndim)]
            p[i] += 1
            I = 0 if i == 0 else np.sum(self._FS_.num_basis._3dCSCG_2Form[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._3dCSCG_2Form[1][i]).reshape(p, order='F') + I,)
        return _ln_
    
    @property
    def _3dCSCG_3Form(self):
        """Return the mesh-element-wise local numbering of 3-form."""
        return [np.arange(self._FS_.num_basis._3dCSCG_3Form[0]).reshape(self._FS_.p, order='F'), ]