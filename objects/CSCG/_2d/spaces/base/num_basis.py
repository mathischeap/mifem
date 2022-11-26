# -*- coding: utf-8 -*-
import numpy as np
from components.freeze.main import FrozenOnly


class NumBasis(FrozenOnly):
    """"""

    def __init__(self, FS):
        """"""
        assert FS.ndim == 2, " <NumBasis> "
        self._FS_ = FS
        self._freeze_self_()



    @property
    def _2dCSCG_0Form_Inner(self):
        """ """
        _basis_ = 1
        for p_i in self._FS_.p:
            _basis_ *= p_i + 1
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_

    @property
    def _2dCSCG_1Form_Inner(self):
        """ """
        _basis_components_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j] + 1 for j in range(self._FS_.ndim)]
            p[i] -= 1
            _basis_components_ += (np.prod(p),)
        _basis_ = np.sum(_basis_components_)
        return _basis_, _basis_components_

    @property
    def _2dCSCG_2Form_Inner(self):
        """ """
        _basis_ = np.prod(self._FS_.p)
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_




    @property
    def _2dCSCG_0Form_Outer(self):
        """ """
        _basis_ = 1
        for p_i in self._FS_.p:
            _basis_ *= p_i + 1
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_

    @property
    def _2dCSCG_1Form_Outer(self):
        """ """
        _basis_components_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j] for j in range(self._FS_.ndim)]
            p[i] += 1
            _basis_components_ += (np.prod(p),)
        _basis_ = np.sum(_basis_components_)
        return _basis_, _basis_components_

    @property
    def _2dCSCG_2Form_Outer(self):
        """ """
        _basis_ = np.prod(self._FS_.p)
        _basis_components_ = (_basis_,)
        return _basis_, _basis_components_




    @property
    def _2dCSCG_0Trace_Inner(self):
        """ """
        p = self._FS_.p
        _basis_ = 2 * ((p[1] + 1) + (p[0] + 1))
        _basis_components_ = {'U': (p[1] + 1,), 'D': (p[1] + 1,),
                              'L': (p[0] + 1,), 'R': (p[0] + 1,)}
        _basis_onsides_ = {'U': p[1] + 1, 'D': p[1] + 1,
                           'L': p[0] + 1, 'R': p[0] + 1}
        return _basis_, _basis_components_, _basis_onsides_

    @property
    def _2dCSCG_1Trace_Inner(self):
        p = self._FS_.p
        _basis_ = 2 * (p[1] + p[0])
        _basis_components_ = {'U': (p[1],), 'D': (p[1],),
                              'L': (p[0],), 'R': (p[0],)}
        _basis_onsides_ = {'U': p[1], 'D': p[1],
                           'L': p[0], 'R': p[0]}
        return _basis_, _basis_components_, _basis_onsides_




    @property
    def _2dCSCG_0Trace_Outer(self):
        """ """
        p = self._FS_.p
        _basis_ = 2 * ((p[1] + 1) + (p[0] + 1))
        _basis_components_ = {'U': (p[1] + 1,), 'D': (p[1] + 1,),
                              'L': (p[0] + 1,), 'R': (p[0] + 1,)}
        _basis_onsides_ = {'U': p[1] + 1, 'D': p[1] + 1,
                           'L': p[0] + 1, 'R': p[0] + 1}
        return _basis_, _basis_components_, _basis_onsides_

    @property
    def _2dCSCG_1Trace_Outer(self):
        p = self._FS_.p
        _basis_ = 2 * (p[1] + p[0])
        _basis_components_ = {'U': (p[1],), 'D': (p[1],),
                              'L': (p[0],), 'R': (p[0],)}
        _basis_onsides_ = {'U': p[1], 'D': p[1],
                           'L': p[0], 'R': p[0]}
        return _basis_, _basis_components_, _basis_onsides_