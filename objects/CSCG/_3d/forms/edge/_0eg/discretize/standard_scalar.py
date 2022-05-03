# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class _3dCSCG_Edge0Form_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()


    def __call__(self, update_cochain=True, target='func'):
        """Discretize the standard _3dCSCG_ScalarField to a 0-edge-form

        'locally full local EEW cochain' means the cochain is a dict whose keys are edge-element
        numbers and values are edge-element-wise local cochains.

        """
        SELF = self._ef_

        nodes = SELF.space.nodes

        local_EEW = dict()
        if target == 'func':
            FUNC = SELF.func.body[0]
        elif target == 'BC':
            FUNC = SELF.BC.body[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_0Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")
        for i in SELF.mesh.edge.elements: # go through all local edge-elements
            element = SELF.mesh.edge.elements[i] # edge-element
            mesh_element = element.CHARACTERISTIC_element
            corner_edge = element.CHARACTERISTIC_corner_edge
            xyz = element.coordinate_transformation.mapping(*nodes,
                                                            from_element=mesh_element,
                                                            corner_edge=corner_edge)

            local_EEW[i] = FUNC(*xyz)
        # isKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: SELF.cochain.local_EEW = local_EEW
        # ...
        return 'locally full local EEW cochain', local_EEW