# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('/')


from components.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.trace._0tr.discretize.vector.standard.flux import _3dCSCG_0Trace_Discretize_StandardVector_Flux
from objects.CSCG._3d.forms.trace._0tr.discretize.scalar.standard import _3dCSCG_0Trace_Discretize_StandardScalar
from objects.CSCG._3d.forms.trace._0tr.discretize.scalar.trace_element_wise import \
    _3dCSCG_0Trace_Discretize_TEW_Scalar

class _3dCSCG_0Trace_Discretize(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._standard_vector_flux_ = _3dCSCG_0Trace_Discretize_StandardVector_Flux(tf)
        self._standard_scalar_ = _3dCSCG_0Trace_Discretize_StandardScalar(tf)
        self._trace_element_wise_scalar_ = _3dCSCG_0Trace_Discretize_TEW_Scalar(tf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param str target:
        :return: The cochain corresponding to the particular discretization scheme.
        """
        SELF = self._tf_
        if target == 'func':
            if SELF.CF.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.CF.ftype == 'standard':
                    return self._standard_scalar_(
                        target = 'func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.CF.ftype}.')

            elif SELF.CF.__class__.__name__ == '_3dCSCG_VectorField':
                if SELF.CF.ftype == 'standard':
                    return self._standard_vector_flux_(
                        target='func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.CF.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target func) '
                                          f'discretize {SELF.CF.body.__class__}.')


        elif target == 'BC':
            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.BC.CF.ftype == 'standard':
                    return self._standard_scalar_(
                        target = 'BC',
                        update_cochain=False)
                elif SELF.BC.CF.ftype == 'trace-element-wise':
                    return self._trace_element_wise_scalar_(
                        target = 'BC',
                        update_cochain=False)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target BC) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.BC.CF.ftype}.')
            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target BC) '
                                          f'discretize {SELF.BC.CF.__class__}.')
        else:
            raise NotImplementedError()





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_0_trace\discretize\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)