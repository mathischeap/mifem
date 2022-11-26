# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('../')

import numpy as np
from components.freeze.base import FrozenOnly



class _3dCSCG_0Trace_Discretize_StandardVector_Flux(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, target='func', update_cochain=True):
        """We will discretize the standard vector field (norm flux) to all trace elements.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
        """
        SELF = self._tf_

        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        nodes = SELF.space.nodes
        nx, ny, nz = nodes
        nodes_NS = np.meshgrid(ny, nz, indexing='ij')
        nodes_WE = np.meshgrid(nx, nz, indexing='ij')
        nodes_BF = np.meshgrid(nx, ny, indexing='ij')

        if target == 'func':
            assert SELF.CF is not None, f"No func.body!"
            fx, fy, fz = SELF.CF.___DO_evaluate_func_at_time___()
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for i in SELF.mesh.trace.elements:
            te = SELF.mesh.trace.elements[i]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                x, y, z = te.coordinate_transformation.mapping(*nodes_NS, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_NS)

            elif ele_side in 'WE':
                x, y, z = te.coordinate_transformation.mapping(*nodes_WE, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_WE)
            elif ele_side in 'BF':
                x, y, z = te.coordinate_transformation.mapping(*nodes_BF, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_BF)
            else:
                raise Exception()

            f0, f1, f2 = fx(x, y, z), fy(x, y, z), fz(x, y, z)

            f = f0 * ouv[0] + f1 * ouv[1] + f2 * ouv[2]

            te_primal_local = f.ravel('F')

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[i] = te_primal_local

        if update_cochain: SELF.cochain.local_TEW = local_TEW

        return 'locally full local TEW cochain', local_TEW





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_0_trace\discretize\vector\standard_flux.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)