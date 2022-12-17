# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/13/2022 4:00 PM
"""
from components.freeze.base import FrozenOnly


class _3dCSCG_0ltf_Discretize_TraceElementWise(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, update_cochain=False, target='BC'):
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
        if target == 'BC':
            assert update_cochain is False, f"Wo do not update cochain when targeting at BC!"
            FUNC = self._ltf_.BC.CF.do.evaluate_func_at_time() # evaluate at current time
        elif target == 'func':
            FUNC = self._ltf_.CF.do.evaluate_func_at_time() # evaluate at current time
        else:
            raise NotImplementedError(f"target={target} is not implemented.")

        s_nodes = self._ltf_.space.nodes

        mesh = self._ltf_.mesh
        trace_elements = mesh.trace.elements
        cochainLocal: dict[int] = dict()

        for e in FUNC:
            te = trace_elements[e]

            if target == 'BC' and not te.whether.on_mesh_boundary:
                # when target at BC, we skip internal trace-elements.
                continue
            else:
                pass

            local_cochain = FUNC[e](*s_nodes)[1][0].ravel('F')

            positions = te.positions

            for pos in positions:

                if pos[0].isnumeric():
                    # we find a mesh-element position for this trace-element
                    element = int(pos[:-1])

                    if element in mesh.elements: # if the mesh-element is local.

                        side = pos[-1] # find the side

                        if element not in cochainLocal:
                            cochainLocal[element]: dict[str] = dict()
                        else:
                            pass

                        cochainLocal[element][side] = local_cochain

                    else:
                        pass

                else:
                    pass

        if update_cochain:
            self._ltf_.cochain.local_ESW = cochainLocal
        else:
            pass

        return 'mesh-element-side-wise local cochain', cochainLocal
