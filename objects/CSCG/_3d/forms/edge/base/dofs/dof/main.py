# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly




class _3dCSCG_Edge_forms_DOF(FrozenOnly):
    """A dof of a edge form."""
    def __init__(self, dofs, i):
        """"""
        # we first check if dof #`i` is a local dof, if not, raise Error.
        ELEMENTS, INDICES = dofs.___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(i)
        self._local_positions_ = list()
        for E, I in zip(ELEMENTS, INDICES):
            self._local_positions_.append((E, I))
        self._i_ = i # I am the #i dof.
        self._dofs_ = dofs
        self._ef_ = dofs._ef_
        self._freeze_self_()