# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:50 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from screws.freeze.main import FrozenOnly

class LocalNumbering(FrozenOnly):
    """ """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <LocalNumbering> "
        self._FS_ = FS
        self._freeze_self_()


    @property
    def _3dCSCG_0Node(self):
        """Return the node-element-wise (NOT mesh-element-wise) local numbering of 0-node-form.

        Sequence: ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF']
        """
        LN = [0,]
        NEW__LN = \
            {'NWB': LN,
             'SWB': LN,
             'NEB': LN,
             'SEB': LN,
             'NWF': LN,
             'SWF': LN,
             'NEF': LN,
             'SEF': LN}
        return NEW__LN

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








    @property
    def _3dCSCG_0Tr(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 0-Tr-form."""
        p = self._FS_.p
        raise NotImplementedError()

    @property
    def _3dCSCG_1Tr(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 1-Tr-form."""
        p = self._FS_.p
        raise NotImplementedError()

    @property
    def _3dCSCG_2Tr(self):
        """Return the trace-element-wise (NOT mesh-element-wise) local numbering of 2-Tr-form."""
        p = self._FS_.p

        num_NS = p[1]*p[2]
        num_WE = p[0]*p[2]
        num_BF = p[0]*p[1]
        TEW__LN = {
            'N': (np.arange(num_NS).reshape((p[1], p[2]), order='F'),),
            'S': (np.arange(num_NS).reshape((p[1], p[2]), order='F'),),
            'W': (np.arange(num_WE).reshape((p[0], p[2]), order='F'),),
            'E': (np.arange(num_WE).reshape((p[0], p[2]), order='F'),),
            'B': (np.arange(num_BF).reshape((p[0], p[1]), order='F'),),
            'F': (np.arange(num_BF).reshape((p[0], p[1]), order='F'),)
            }

        # the gathering matrix to gather local Trace-Element-Wise to Mesh-Element-Wise
        local_gathering_matrix = {
            'N': np.arange(num_NS),
            'S': np.arange(num_NS) + num_NS,
            'W': np.arange(num_NS) + num_NS * 2,
            'E': np.arange(num_NS) + num_NS * 2 + num_WE,
            'B': np.arange(num_NS) + num_NS * 2 + num_WE * 2,
            'F': np.arange(num_NS) + num_NS * 2 + num_WE * 2 + num_BF,
        }

        return TEW__LN, local_gathering_matrix










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