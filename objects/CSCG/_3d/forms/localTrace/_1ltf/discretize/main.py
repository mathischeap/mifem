# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:24 PM
"""
from components.freeze.main import FrozenOnly

from objects.CSCG._3d.forms.localTrace._1ltf.discretize.vector.trace_element_wise import _3dCSCG_1LocalTrace_Discretize_TraceElementWiseVector
from objects.CSCG._3d.forms.localTrace._1ltf.discretize.vector.standard.T_perp import _3dCSCG_1LocalTrace_Discretize_StandardVector_T_perp
from objects.CSCG._3d.forms.localTrace._1ltf.discretize.vector.standard.T_para import _3dCSCG_1LocalTrace_Discretize_StandardVector_T_para


class _3dCSCG_1LocalTrace_Discretize(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._standard_vector_T_para_ = _3dCSCG_1LocalTrace_Discretize_StandardVector_T_para(ltf)
        self._standard_vector_T_perp_ = _3dCSCG_1LocalTrace_Discretize_StandardVector_T_perp(ltf)
        self._trace_element_wise_vector_ = _3dCSCG_1LocalTrace_Discretize_TraceElementWiseVector(ltf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', component='T_para'):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param target:
        :param component:
        :return: The cochain corresponding to the particular discretization scheme.
        """
        SELF = self._ltf_


        if target == 'func':

            if SELF.CF.__class__.__name__ == '_3dCSCG_VectorField':


                if SELF.CF.ftype == 'standard':

                    if component == 'T_para':

                        return self._standard_vector_T_para_(
                            update_cochain=update_cochain, target='func')

                    elif component == 'T_perp':

                        return self._standard_vector_T_perp_(
                            update_cochain=update_cochain, target='func')

                    else:
                        raise Exception(f"1-ltf discretization targeting vector func of standard type. "
                                        f"I cannot discretize component = {component}. It should be either"
                                        f"'T_para' (parallel trace) or 'T_perp' (perpendicular trace).")


                elif SELF.CF.ftype == 'trace-element-wise':
                    # we do not care this trace-element-wise vector is T_para or T_perp vector, we just discretize it to the trace.
                    return self._trace_element_wise_vector_(
                        update_cochain=update_cochain, target='func')


                else:
                    raise Exception(f'3dCSCG 1-ltf can not (target func) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.CF.ftype}.')



            else:
                raise NotImplementedError(f'3dCSCG 1-ltf can not (target func) '
                                          f'discretize {SELF.CF.__class__}.')

        elif target == 'BC':

            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.BC.CF.ftype == 'trace-element-wise':
                    # we do not care this trace-element-wise vector is T_para or T_perp vector, we just discretize it to the trace.
                    return self._trace_element_wise_vector_(
                        update_cochain=False, target='BC')

                else:
                    raise Exception(f'3dCSCG 1-lft can not (target BC) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.BC.CF.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 1-ltf can not (target BC) '
                                          f'discretize {SELF.BC.CF__class__}.')


        else:
            raise NotImplementedError(f"target={target} not implemented "
                                      f"for 3d CSCG 1-ltf form discretization.")
