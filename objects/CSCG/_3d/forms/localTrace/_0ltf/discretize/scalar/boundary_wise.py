# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np

class _3dCSCG_0ltf_Discretize_BoundaryWise(FrozenOnly):
    def __init__(self, ltf):
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """

        'Boundary only local cochain' means we return a dict, its keys are mesh-element numbers,
        its values are also dictionaries whose keys are mesh-element-side names, like 'N', 'S' and
        so on, and values are the mesh-element-side(trace-element)-wise local cochains. For example
        cochainLocal = {
                1: {'N': [4, 3, 1, 1.5, ...], 'W': [...]},
                23: {...},
                ...}
        We know we have cochains for mesh-element #1, #23, ..., and for mesh element #1, we have
        local cochain on its North side and West side.
        """
        if target == 'func':
            FUNC = self._ltf_.CF.do.evaluate_func_at_time()
        elif target == 'BC':
            FUNC = self._ltf_.BC.CF.do.evaluate_func_at_time()
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise Exception()

        Nsx, Nsy, Nsz = self._ltf_.space.nodes
        nodes_NS = [np.meshgrid(Nsy, Nsz, indexing='ij')[i].ravel('F') for i in range(2)]
        nodes_WE = [np.meshgrid(Nsx, Nsz, indexing='ij')[i].ravel('F') for i in range(2)]
        nodes_BF = [np.meshgrid(Nsx, Nsy, indexing='ij')[i].ravel('F') for i in range(2)]

        ND = {
            'N': nodes_NS,
            'S': nodes_NS,
            'W': nodes_WE,
            'E': nodes_WE,
            'B': nodes_BF,
            'F': nodes_BF
        }

        RANGE_element_sides = self._ltf_.mesh.boundaries.range_of_element_sides
        cochainLocal = dict()
        Tmap = self._ltf_.mesh.trace.elements.map
        NSWEBF = 'NSWEBF'

        for bn in FUNC:
            func_bn = FUNC[bn][0]
            element_sides = RANGE_element_sides[bn]
            elements, sides = list(), list()
            for element_side in element_sides:
                element = int(element_side[:-1])
                side = element_side[-1]
                elements.append(element)
                sides.append(side)
            for i, side in zip(elements, sides):
                te = Tmap[i][NSWEBF.index(side)]
                xyz = self._ltf_.mesh.trace.elements[te].coordinate_transformation.mapping(
                    *ND[side], from_element=i, side=side, parse_3_1d_eps=False)
                local_cochain = func_bn(*xyz)

                if i not in cochainLocal:
                    cochainLocal[i] = dict()
                else:
                    pass

                cochainLocal[i][side] = local_cochain

        if update_cochain:
            self._ltf_.cochain.local_ESW = cochainLocal
        else:
            pass

        return 'mesh-element-side-wise local cochain', cochainLocal