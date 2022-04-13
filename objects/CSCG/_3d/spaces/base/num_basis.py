# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:49 2019
Aerodynamics, AE
TU Delft
"""


import numpy as np
from screws.freeze.main import FrozenOnly




class NumBasis(FrozenOnly):
    """ """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <NumBasis> "
        self._FS_ = FS
        self._freeze_self_()


    @property
    def _3dCSCG_0Node(self):
        """ """
        _basis_ = 8
        _basis_components_ = {'NWB': (1,),
                              'SWB': (1,),
                              'NEB': (1,),
                              'SEB': (1,),
                              'NWF': (1,),
                              'SWF': (1,),
                              'NEF': (1,),
                              'SEF': (1,)}
        return _basis_, _basis_components_




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
        _basis_ = 2 * ((p[1]+1)*(p[2]+1)+(p[0]+1)*(p[2]+1)+(p[0]+1)*(p[1]+1))
        _basis_components_ = {'N': ((p[1]+1)*(p[2]+1),), 'S': ((p[1]+1)*(p[2]+1),),
                              'W': ((p[0]+1)*(p[2]+1),), 'E': ((p[0]+1)*(p[2]+1),),
                              'B': ((p[0]+1)*(p[1]+1),), 'F': ((p[0]+1)*(p[1]+1),)}
        _basis_onsides_ = {'N': (p[1]+1)*(p[2]+1), 'S': (p[1]+1)*(p[2]+1),
                           'W': (p[0]+1)*(p[2]+1), 'E': (p[0]+1)*(p[2]+1),
                           'B': (p[0]+1)*(p[1]+1), 'F': (p[0]+1)*(p[1]+1)}
        return _basis_, _basis_components_, _basis_onsides_
    
    @property
    def _3dCSCG_1Trace(self):
        """ """
        p = self._FS_.p
        _basis_ = 2 * (p[1]*(p[2]+1) + (p[1]+1)*p[2] + 
                       p[0]*(p[2]+1) + (p[0]+1)*p[2] + 
                       p[0]*(p[1]+1) + (p[0]+1)*p[1])
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
        return _basis_, _basis_components_, _basis_onsides_








    @property
    def _3dCSCG_0Tr(self):
        p = self._FS_.p
        _basis_ = (p[0] + 1) * (p[1] + 1) * (p[2] + 1) - (p[0]-1) * (p[1]-1) * (p[2]-1)

        _basis_components_ = {'N': ((p[1]+1)*(p[2]+1),), 'S': ((p[1]+1)*(p[2]+1),),
                              'W': ((p[0]+1)*(p[2]+1),), 'E': ((p[0]+1)*(p[2]+1),),
                              'B': ((p[0]+1)*(p[1]+1),), 'F': ((p[0]+1)*(p[1]+1),)}

        _basis_onsides_ = {'N': (p[1]+1)*(p[2]+1), 'S': (p[1]+1)*(p[2]+1),
                           'W': (p[0]+1)*(p[2]+1), 'E': (p[0]+1)*(p[2]+1),
                           'B': (p[0]+1)*(p[1]+1), 'F': (p[0]+1)*(p[1]+1)}
        return _basis_, _basis_components_, _basis_onsides_

    @property
    def _3dCSCG_1Tr(self):
        """ """
        p = self._FS_.p
        _basis_ = p[0] * (p[1] + 1) * (p[2] + 1) - p[0] * (p[1] - 1) * (p[2] - 1) + \
                  (p[0] + 1) * p[1] * (p[2] + 1) - (p[0] - 1) * p[1] * (p[2] - 1) + \
                  (p[0] + 1) * (p[1] + 1) * p[2] - (p[0] - 1) * (p[1] - 1) * p[2]

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

        return _basis_, _basis_components_, _basis_onsides_

    @property
    def _3dCSCG_2Tr(self):
        p = self._FS_.p

        _basis_ = 2 * (p[1]*p[2] + p[0]*p[2] + p[0]*p[1])

        _basis_components_ = {'N': (p[1]*p[2],), 'S': (p[1]*p[2],),
                              'W': (p[0]*p[2],), 'E': (p[0]*p[2],),
                              'B': (p[0]*p[1],), 'F': (p[0]*p[1],)}
        _basis_onsides_ = {'N': p[1]*p[2], 'S': p[1]*p[2],
                           'W': p[0]*p[2], 'E': p[0]*p[2],
                           'B': p[0]*p[1], 'F': p[0]*p[1]}

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