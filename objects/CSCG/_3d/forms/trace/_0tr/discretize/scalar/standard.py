# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('/')
import numpy as np

from screws.freeze.base import FrozenOnly



class _3dCSCG_0Trace_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, target='func', update_cochain=True):
        """We will discretize the standard scalar field to all trace elements.

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
            _lf_ = SELF.CF.___DO_evaluate_func_at_time___()[0]
        elif target == 'BC':
            assert SELF.BC.CF is not None, f"No BC.body!"
            _lf_ = SELF.BC.CF.___DO_evaluate_func_at_time___()[0]
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for i in SELF.mesh.trace.elements:
            te = SELF.mesh.trace.elements[i]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                x, y, z = te.coordinate_transformation.mapping(*nodes_NS, from_element=ele, side=ele_side)
            elif ele_side in 'WE':
                x, y, z = te.coordinate_transformation.mapping(*nodes_WE, from_element=ele, side=ele_side)
            elif ele_side in 'BF':
                x, y, z = te.coordinate_transformation.mapping(*nodes_BF, from_element=ele, side=ele_side)
            else:
                raise Exception()

            te_primal_local = _lf_(x, y, z).ravel('F')

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[i] = te_primal_local

        if update_cochain: SELF.cochain.local_TEW = local_TEW

        return 'locally full local TEW cochain', local_TEW





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_0_trace\discretize\scalar\standard.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)