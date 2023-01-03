# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:24 PM
"""
from components.freeze.main import FrozenOnly

from objects.CSCG._3d.forms.localTrace._0ltf.discretize.scalar.boundary_wise import _3dCSCG_0ltf_Discretize_BoundaryWise
from objects.CSCG._3d.forms.localTrace._0ltf.discretize.scalar.standard import _3dCSCG_0ltf_Discretize_Standard
from objects.CSCG._3d.forms.localTrace._0ltf.discretize.scalar.trace_element_wise import \
    _3dCSCG_0ltf_Discretize_TraceElementWise


class _3dCSCG_0LocalTrace_Discretize(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._standard_ = _3dCSCG_0ltf_Discretize_Standard(ltf)
        self._boundary_wise_ = _3dCSCG_0ltf_Discretize_BoundaryWise(ltf)
        self._trace_element_wise_ = _3dCSCG_0ltf_Discretize_TraceElementWise(ltf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.

        """
        SELF = self._ltf_
        if target == 'func':
            if SELF.CF.__class__.__name__ == '_3dCSCG_ScalarField':

                if SELF.CF.ftype == 'standard':
                    return self._standard_(update_cochain=update_cochain, target='func')

                elif SELF.CF.ftype == 'boundary-wise':
                    return self._boundary_wise_(target='func')

                else:
                    raise NotImplementedError(f"3dCSCG 0-ltf cannot (target func) discretize "
                                              f"_3dCSCG_ScalarField of ftype={SELF.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-ltf can not (target func) discretize '
                                          f'{SELF.CF.__class__}.')

        elif target == 'BC':

            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_ScalarField':

                if SELF.BC.CF.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self._standard_(update_cochain=False, target='BC')

                elif SELF.BC.CF.ftype == "boundary-wise":
                    # we will always not update cochain & and always set target to be "BC"
                    return self._boundary_wise_(target='BC')

                elif SELF.BC.CF.ftype == "trace-element-wise":
                    return self._trace_element_wise_(update_cochain=False, target='BC')

                else:
                    raise NotImplementedError(f"3dCSCG 0-ltf cannot (target BC) discretize "
                                              f"_3dCSCG_ScalarField of ftype={SELF.BC.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-ltf can not (target BC) '
                                          f'discretize {SELF.BC.CF.__class__}.')

        else:
            raise NotImplementedError(f"3dCSCG 0-ltf cannot discretize while targeting at {target}.")
