# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np


class _2dCSCG_S0F_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        nodes = list(np.meshgrid(*self._sf_.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(2)]
        cochainLocal = dict()

        if target == 'func':
            FUNC = self._sf_.CF.___DO_evaluate_func_at_time___()[0]
        elif target == 'BC':
            FUNC = self._sf_.BC.CF.___DO_evaluate_func_at_time___()[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_2d_outer_0Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")

        for i in self._sf_.mesh.elements:
            element = self._sf_.mesh.elements[i]
            xyz = element.coordinate_transformation.mapping(*nodes)
            cochainLocal[i] = FUNC(*xyz)
        # isKronecker? ...
        if not self._sf_.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self._sf_.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal