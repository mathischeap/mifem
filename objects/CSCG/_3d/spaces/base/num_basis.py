# -*- coding: utf-8 -*-
"""INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:49 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from components.freeze.main import FrozenOnly


class NumBasis(FrozenOnly):
    """ """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <NumBasis> "
        self._FS_ = FS
        self._freeze_self_()

    @property
    def _3dCSCG_0Edge(self):
        """ """
        p = self._FS_.p
        _basis_ = 4 * (p[0] + 1 + p[1] + 1 + p[2] + 1)
        _basis_components_ = {'WB': (p[0] + 1,),
                              'EB': (p[0] + 1,),
                              'WF': (p[0] + 1,),
                              'EF': (p[0] + 1,),
                              'NB': (p[1] + 1,),
                              'SB': (p[1] + 1,),
                              'NF': (p[1] + 1,),
                              'SF': (p[1] + 1,),
                              'NW': (p[2] + 1,),
                              'SW': (p[2] + 1,),
                              'NE': (p[2] + 1,),
                              'SE': (p[2] + 1,)}
        return _basis_, _basis_components_

    @property
    def _3dCSCG_1Edge(self):
        """ """
        p = self._FS_.p
        _basis_ = 4 * (p[0] + p[1] + p[2])
        _basis_components_ = {'WB': (p[0],),
                              'EB': (p[0],),
                              'WF': (p[0],),
                              'EF': (p[0],),
                              'NB': (p[1],),
                              'SB': (p[1],),
                              'NF': (p[1],),
                              'SF': (p[1],),
                              'NW': (p[2],),
                              'SW': (p[2],),
                              'NE': (p[2],),
                              'SE': (p[2],)}
        return _basis_, _basis_components_

    @property
    def _3dCSCG_0Trace(self):
        """ """
        p = self._FS_.p
        _basis_hybrid_True_ = 2 * ((p[1]+1)*(p[2]+1)+(p[0]+1)*(p[2]+1)+(p[0]+1)*(p[1]+1))

        px, py, pz = p
        _basis_hybrid_False_ = \
            _basis_hybrid_True_ \
            - 8 * 2 \
            - (px - 1) * 4 \
            - (py - 1) * 4 \
            - (pz - 1) * 4

        _basis_components_ = {'N': ((p[1]+1)*(p[2]+1),), 'S': ((p[1]+1)*(p[2]+1),),
                              'W': ((p[0]+1)*(p[2]+1),), 'E': ((p[0]+1)*(p[2]+1),),
                              'B': ((p[0]+1)*(p[1]+1),), 'F': ((p[0]+1)*(p[1]+1),)}
        _basis_onsides_ = {'N': (p[1]+1)*(p[2]+1), 'S': (p[1]+1)*(p[2]+1),
                           'W': (p[0]+1)*(p[2]+1), 'E': (p[0]+1)*(p[2]+1),
                           'B': (p[0]+1)*(p[1]+1), 'F': (p[0]+1)*(p[1]+1)}

        _basis_ = {
            True: _basis_hybrid_True_,
            False: _basis_hybrid_False_,
        }

        return _basis_, _basis_components_, _basis_onsides_
    
    @property
    def _3dCSCG_1Trace(self):
        """ """
        p = self._FS_.p
        _basis_hybrid_True_ = 2 * (
            p[1]*(p[2]+1) + (p[1]+1)*p[2] +
            p[0]*(p[2]+1) + (p[0]+1)*p[2] +
            p[0]*(p[1]+1) + (p[0]+1)*p[1]
        )

        px, py, pz = p
        _basis_hybrid_False_ = \
            _basis_hybrid_True_ \
            - px * 4 \
            - py * 4 \
            - pz * 4

        _basis_components_ = {'N': (p[1]*(p[2]+1), (p[1]+1)*p[2]), 
                              'S': (p[1]*(p[2]+1), (p[1]+1)*p[2]), 
                              'W': (p[0]*(p[2]+1), (p[0]+1)*p[2]), 
                              'E': (p[0]*(p[2]+1), (p[0]+1)*p[2]), 
                              'B': (p[0]*(p[1]+1), (p[0]+1)*p[1]), 
                              'F': (p[0]*(p[1]+1), (p[0]+1)*p[1])}
        _basis_onsides_ = {'N': p[1]*(p[2]+1) + (p[1]+1)*p[2], 
                           'S': p[1]*(p[2]+1) + (p[1]+1)*p[2], 
                           'W': p[0]*(p[2]+1) + (p[0]+1)*p[2],
                           'E': p[0]*(p[2]+1) + (p[0]+1)*p[2],
                           'B': p[0]*(p[1]+1) + (p[0]+1)*p[1],
                           'F': p[0]*(p[1]+1) + (p[0]+1)*p[1]}

        _basis_ = {
            True: _basis_hybrid_True_,
            False: _basis_hybrid_False_,
        }

        return _basis_, _basis_components_, _basis_onsides_
    
    @property
    def _3dCSCG_2Trace(self):
        p = self._FS_.p
        _basis_ = 2 * (p[1]*p[2] + p[0]*p[2] + p[0]*p[1])
        _basis_components_ = {'N': (p[1]*p[2],), 'S': (p[1]*p[2],),
                              'W': (p[0]*p[2],), 'E': (p[0]*p[2],),
                              'B': (p[0]*p[1],), 'F': (p[0]*p[1],)}
        _basis_onsides_ = {'N': p[1]*p[2], 'S': p[1]*p[2],
                           'W': p[0]*p[2], 'E': p[0]*p[2],
                           'B': p[0]*p[1], 'F': p[0]*p[1]}
        _basis_ = {
            True: _basis_,
            False: _basis_,
        }
        return _basis_, _basis_components_, _basis_onsides_

    @property
    def _3dCSCG_0Form(self):
        """ """
        _basis_ = 1
        for p_i in self._FS_.p:
            _basis_ *= p_i + 1
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_
    
    @property
    def _3dCSCG_1Form(self):
        """ """
        _basis_components_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[i]+1 for i in range(self._FS_.ndim)]
            p[i] -= 1
            _basis_components_ += (np.prod(p),)
        _basis_ = np.sum(_basis_components_)
        return _basis_, _basis_components_
    
    @property
    def _3dCSCG_2Form(self):
        """ """
        _basis_components_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[i] for i in range(self._FS_.ndim)]
            p[i] += 1
            _basis_components_ += (np.prod(p),)
        _basis_ = np.sum(_basis_components_)
        return _basis_, _basis_components_
    
    @property
    def _3dCSCG_3Form(self):
        """ """
        _basis_ = np.prod(self._FS_.p)
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_

    @property
    def _3dCSCG_0LocalTrace(self):
        px, py, pz = self._FS_.p
        px += 1
        py += 1
        pz += 1
        basis = {
            True: 2 * (px * py + py * pz + pz * px),
            False: self._3dCSCG_0Trace[0][False],
        }
        return basis, {
            'N': py * pz,
            'S': py * pz,
            'W': px * pz,
            'E': px * pz,
            'B': px * py,
            'F': px * py
        }

    @property
    def _3dCSCG_2LocalTrace(self):
        px, py, pz = self._FS_.p
        basis = {
            True: 2 * (px * py + py * pz + pz * px),
            False: self._3dCSCG_2Trace[0][False],
        }
        return basis, {
            'N': py * pz,
            'S': py * pz,
            'W': px * pz,
            'E': px * pz,
            'B': px * py,
            'F': px * py
        }

    @property
    def _3dCSCG_1LocalTrace(self):
        px, py, pz = self._FS_.p
        num_x_face = (py + 1) * pz + (pz + 1) * py
        num_y_face = (px + 1) * pz + (pz + 1) * px
        num_z_face = (px + 1) * py + (py + 1) * px
        basis = {
            True: 2 * (num_x_face + num_y_face + num_z_face),
            False: self._3dCSCG_1Trace[0][False],
        }
        return basis, {
            'N': num_x_face,
            'S': num_x_face,
            'W': num_y_face,
            'E': num_y_face,
            'B': num_z_face,
            'F': num_z_face,
        }
