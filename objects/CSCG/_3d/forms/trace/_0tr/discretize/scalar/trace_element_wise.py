# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class _3dCSCG_0Trace_Discretize_TEW_Scalar(FrozenOnly):
    """"""
    def __init__(self, tf):
        """"""
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, target='func', update_cochain=True):
        """"""
        SELF = self._tf_
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        nodes = SELF.space.nodes
        nx, ny, nz = nodes

        if target == 'func':
            assert SELF.CF is not None, f"No func.body!"
            TEW_FUNC = SELF.CF.___DO_evaluate_func_at_time___()
        elif target == 'BC':
            assert SELF.BC.CF is not None, f"No BC.body!"
            TEW_FUNC = SELF.BC.CF.___DO_evaluate_func_at_time___()
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for i in SELF.mesh.trace.elements:

            te_primal_local = TEW_FUNC[i](nx, ny, nz)[1][0].ravel('F')

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[i] = te_primal_local

        if update_cochain: SELF.cochain.local_TEW = local_TEW

        return 'locally full local TEW cochain', local_TEW