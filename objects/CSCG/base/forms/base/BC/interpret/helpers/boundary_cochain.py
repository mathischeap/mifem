# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/15/2022 4:42 PM
"""
from scipy.sparse import csr_matrix
from components.freeze.main import FrozenOnly


class CSCG_FORM_BC_Interpret_BoundaryCochain(FrozenOnly):
    """It will follow the s.BC.CF and s.BC.boundaries in real time."""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        self._elements_ = f.mesh.elements
        self.___empty___ = csr_matrix((f.num.basis, 1))
        self._freeze_self_()

    def __call__(self, i):
        """Return the boundary local cochains in real time (following the current BC.CF and BC.boundaries).
        for mesh-element i. For those dofs not locating on the mesh boundary, we return 0 for them.

        Parameters
        ----------
        i

        Returns
        -------

        """
        element = self._elements_[i]
        if element.whether.internal:
            return self.___empty___
        else:
            raise NotImplementedError()