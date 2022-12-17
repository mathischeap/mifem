# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np


class _3dCSCG_Discretize_Standard(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        SELF = self._sf_

        nodes = list(np.meshgrid(*SELF.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(3)]
        cochainLocal = dict()
        if target == 'func':
            FUNC = SELF.CF.do.evaluate_func_at_time()[0]
        elif target == 'BC':
            FUNC = SELF.BC.CF.do.evaluate_func_at_time()[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_0Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")
        for i in SELF.mesh.elements:
            element = SELF.mesh.elements[i]
            xyz = element.coordinate_transformation.mapping(*nodes)
            cochainLocal[i] = FUNC(*xyz)
        # isKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal
