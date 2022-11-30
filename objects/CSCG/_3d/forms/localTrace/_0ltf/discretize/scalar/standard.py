# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np


class _3dCSCG_0ltf_Discretize_Standard(FrozenOnly):
    def __init__(self, ltf):
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        if target == 'func':
            FUNC = self._ltf_.CF.do.evaluate_func_at_time()[0] # evaluate at current time
        elif target == 'BC':
            FUNC = self._ltf_.BC.CF.do.evaluate_func_at_time()[0] # evaluate at current time
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"0ltf.discretize_standard does not work for target={target}.")

        cochainLocal = dict()

        Nsx, Nsy, Nsz = self._ltf_.space.nodes
        nodes_NS = [np.meshgrid(Nsy, Nsz, indexing='ij')[i].ravel('F') for i in range(2)]
        nodes_WE = [np.meshgrid(Nsx, Nsz, indexing='ij')[i].ravel('F') for i in range(2)]
        nodes_BF = [np.meshgrid(Nsx, Nsy, indexing='ij')[i].ravel('F') for i in range(2)]

        Tmap = self._ltf_.mesh.trace.elements.map
        for i in self._ltf_.mesh.elements:
            tes = Tmap[i]
            xyz_N = self._ltf_.mesh.trace.elements[tes[0]].coordinate_transformation.mapping(
                *nodes_NS, from_element=i, side='N', parse_3_1d_eps=False)
            xyz_S = self._ltf_.mesh.trace.elements[tes[1]].coordinate_transformation.mapping(
                *nodes_NS, from_element=i, side='S', parse_3_1d_eps=False)
            xyz_W = self._ltf_.mesh.trace.elements[tes[2]].coordinate_transformation.mapping(
                *nodes_WE, from_element=i, side='W', parse_3_1d_eps=False)
            xyz_E = self._ltf_.mesh.trace.elements[tes[3]].coordinate_transformation.mapping(
                *nodes_WE, from_element=i, side='E', parse_3_1d_eps=False)
            xyz_B = self._ltf_.mesh.trace.elements[tes[4]].coordinate_transformation.mapping(
                *nodes_BF, from_element=i, side='B', parse_3_1d_eps=False)
            xyz_F = self._ltf_.mesh.trace.elements[tes[5]].coordinate_transformation.mapping(
                *nodes_BF, from_element=i, side='F', parse_3_1d_eps=False)

            cochainLocal[i] = np.concatenate([
                FUNC(*xyz_N),
                FUNC(*xyz_S),
                FUNC(*xyz_W),
                FUNC(*xyz_E),
                FUNC(*xyz_B),
                FUNC(*xyz_F)
            ])

        # isKronecker? ...
        if not self._ltf_.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self._ltf_.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal