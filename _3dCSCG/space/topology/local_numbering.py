# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:50 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from SCREWS.frozen import FrozenOnly

class LocalNumbering(FrozenOnly):
    """ """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <LocalNumbering> "
        self._FS_ = FS
        self._freeze_self_()

    @property
    def _0Trace(self):
        """ """
        p = self._FS_.p
        _local_ = {'N': (np.arange((p[1]+1)*(p[2]+1)).reshape((p[1]+1, p[2]+1), order='F'),),
                   'S': (np.arange((p[1]+1)*(p[2]+1)).reshape((p[1]+1, p[2]+1), order='F'),),
                   'W': (np.arange((p[0]+1)*(p[2]+1)).reshape((p[0]+1, p[2]+1), order='F'),),
                   'E': (np.arange((p[0]+1)*(p[2]+1)).reshape((p[0]+1, p[2]+1), order='F'),),
                   'B': (np.arange((p[0]+1)*(p[1]+1)).reshape((p[0]+1, p[1]+1), order='F'),),
                   'F': (np.arange((p[0]+1)*(p[1]+1)).reshape((p[0]+1, p[1]+1), order='F'),)}
        return _local_
    
    @property
    def _1Trace(self):
        """ """
        p = self._FS_.p
        _local_ = dict()
        _local_['N'] = (np.arange(p[1]*(p[2]+1)).reshape((p[1], p[2]+1), order='F'),
                        np.arange((p[1]+1)*p[2]).reshape((p[1]+1, p[2]), order='F') + p[1]*(p[2]+1))
        _local_['S'] = (np.arange(p[1]*(p[2]+1)).reshape((p[1], p[2]+1), order='F'),
                        np.arange((p[1]+1)*p[2]).reshape((p[1]+1, p[2]), order='F') + p[1]*(p[2]+1))
        _local_['W'] = (np.arange(p[0]*(p[2]+1)).reshape((p[0], p[2]+1), order='F'),
                        np.arange((p[0]+1)*p[2]).reshape((p[0]+1, p[2]), order='F') + p[0]*(p[2]+1))
        _local_['E'] = (np.arange(p[0]*(p[2]+1)).reshape((p[0], p[2]+1), order='F'),
                        np.arange((p[0]+1)*p[2]).reshape((p[0]+1, p[2]), order='F') + p[0]*(p[2]+1))
        _local_['B'] = (np.arange(p[0]*(p[1]+1)).reshape((p[0], p[1]+1), order='F'),
                        np.arange((p[0]+1)*p[1]).reshape((p[0]+1, p[1]), order='F') + p[0]*(p[1]+1))
        _local_['F'] = (np.arange(p[0]*(p[1]+1)).reshape((p[0], p[1]+1), order='F'),
                        np.arange((p[0]+1)*p[1]).reshape((p[0]+1, p[1]), order='F') + p[0]*(p[1]+1))
        return _local_
    
    @property
    def _2Trace(self):
        p = self._FS_.p
        _local_ = {'N': (np.arange(p[1]*p[2]).reshape((p[1], p[2]), order='F'),),
                   'S': (np.arange(p[1]*p[2]).reshape((p[1], p[2]), order='F'),),
                   'W': (np.arange(p[0]*p[2]).reshape((p[0], p[2]), order='F'),),
                   'E': (np.arange(p[0]*p[2]).reshape((p[0], p[2]), order='F'),),
                   'B': (np.arange(p[0]*p[1]).reshape((p[0], p[1]), order='F'),),
                   'F': (np.arange(p[0]*p[1]).reshape((p[0], p[1]), order='F'),)}
        return _local_
    
    @property
    def _0Form(self):
        """ """
        p = [self._FS_.p[i]+1 for i in range(self._FS_.ndim)]
        _ln_ = (np.arange(self._FS_.num_basis._0Form[0]).reshape(*p, order='F'),)
        return _ln_
    
    @property
    def _1Form(self):
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j]+1 for j in range(self._FS_.ndim)]
            p[i] -= 1
            I = 0 if i ==0 else np.sum(self._FS_.num_basis._1Form[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._1Form[1][i]).reshape(p, order='F') + I,)
        return _ln_
    
    @property
    def _2Form(self):
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j] for j in range(self._FS_.ndim)]
            p[i] += 1
            I = 0 if i == 0 else np.sum(self._FS_.num_basis._2Form[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._2Form[1][i]).reshape(p, order='F') + I,)
        return _ln_
    
    @property
    def _3Form(self):
        return [np.arange(self._FS_.num_basis._3Form[0]).reshape(self._FS_.p, order='F'),]
    