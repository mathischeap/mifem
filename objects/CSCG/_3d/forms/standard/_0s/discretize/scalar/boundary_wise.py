# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np

class _3dCSCG_Discretize_BoundaryWise(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()


    def __call__(self):
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
        SELF = self._sf_

        nodes = list(np.meshgrid(*SELF.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(3)]
        FUNC = SELF.BC.CF.do.evaluate_func_at_time()
        RANGE_element_sides = SELF.mesh.boundaries.range_of_element_sides
        cochainLocal = dict()
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
                element = SELF.mesh.elements[i]
                xyz = element.coordinate_transformation.mapping(*nodes)
                local_cochain = func_bn(*xyz)
                local_dofs = SELF.numbering.do.find.local_dofs_on_element_side(side)
                if i not in cochainLocal:
                    cochainLocal[i] = dict()

                cochainLocal[i][side] = local_cochain[local_dofs]

        return 'Boundary only local cochain', cochainLocal