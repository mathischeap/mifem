# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:24 PM
"""
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.localTrace._2ltf.discretize.scalar.standard import _3dCSCG_2ltf_Discretize_Standard



class _3dCSCG_2LocalTrace_Discretize(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._standard_ = _3dCSCG_2ltf_Discretize_Standard(ltf)
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
                else:
                    raise NotImplementedError(f"3dCSCG 2-ltf cannot (target func) discretize "
                                              f"_3dCSCG_ScalarField of ftype={SELF.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-ltf can not (target func) discretize '
                                          f'{SELF.CF.__class__}.')

        elif target == 'BC':
            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.BC.CF.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self._standard_(update_cochain=False, target='BC')

                else:
                    raise NotImplementedError(f"3dCSCG 2-ltf cannot (target BC) discretize "
                                              f"_3dCSCG_ScalarField of ftype={SELF.BC.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 2-ltf can not (target BC) '
                                          f'discretize {SELF.BC.CF.__class__}.')
        else:
            raise NotImplementedError(f"3dCSCG 2-ltf cannot discretize while targeting at {target}.")

