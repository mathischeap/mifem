# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.trace._1tr.discretize.vector.standard.T_perp import _3dCSCG_1Trace_Discretize_StandardVector_T_perp
from objects.CSCG._3d.forms.trace._1tr.discretize.vector.standard.T_para import _3dCSCG_1Trace_Discretize_StandardVector_T_para
from objects.CSCG._3d.forms.trace._1tr.discretize.vector.trace_element_wise import _3dCSCG_1Trace_Discretize_TraceElementWiseVector



class _3dCSCG_1Trace_Discretize(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._standard_vector_T_para_ = _3dCSCG_1Trace_Discretize_StandardVector_T_para(tf)
        self._standard_vector_T_perp_ = _3dCSCG_1Trace_Discretize_StandardVector_T_perp(tf)
        self._trace_element_wise_vector_ = _3dCSCG_1Trace_Discretize_TraceElementWiseVector(tf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', component='T_para', **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param target:
        :param component:
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """
        SELF = self._tf_


        if target == 'func':

            if SELF.CF.__class__.__name__ == '_3dCSCG_VectorField':


                if SELF.CF.ftype == 'standard':

                    if component == 'T_para':

                        return self._standard_vector_T_para_(
                            update_cochain=update_cochain, target='func', **kwargs)

                    elif component == 'T_perp':

                        return self._standard_vector_T_perp_(
                            update_cochain=update_cochain, target='func', **kwargs)

                    else:
                        raise Exception(f"1-trace-form discretization targeting vector func of standard type. "
                                        f"I cannot discretize component = {component}. It should be either"
                                        f"'T_para' (parallel trace) or 'T_perp' (perpendicular trace).")


                elif SELF.CF.ftype == 'trace-element-wise':
                    # we do not care this trace-element-wise vector is T_para or T_perp vector, we just discretize it to the trace.
                    return self._trace_element_wise_vector_(
                        update_cochain=update_cochain, target='func', **kwargs)


                else:
                    raise Exception(f'3dCSCG 1-trace can not (target func) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.CF.ftype}.')



            else:
                raise NotImplementedError(f'3dCSCG 1-trace can not (target func) '
                                          f'discretize {SELF.CF.__class__}.')

        elif target == 'BC':

            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.BC.CF.ftype == 'trace-element-wise':
                    # we do not care this trace-element-wise vector is T_para or T_perp vector, we just discretize it to the trace.
                    return self._trace_element_wise_vector_(
                        update_cochain=False, target='BC', **kwargs)

                else:
                    raise Exception(f'3dCSCG 1-trace can not (target BC) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.BC.CF.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 1-trace can not (target BC) '
                                          f'discretize {SELF.BC.CF__class__}.')


        else:
            raise NotImplementedError(f"target={target} not implemented "
                                      f"for 3d CSCG 1-trace form discretization.")




if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_1_trace\discretize\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)