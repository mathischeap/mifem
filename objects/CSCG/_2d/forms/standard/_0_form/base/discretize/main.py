# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard._0_form.base.discretize.scalar.standard import _2dCSCG_S0F_Discretize_StandardScalar

class _2dCSCG_S0F_Discretize(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._standard_scalar_ = _2dCSCG_S0F_Discretize_StandardScalar(sf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':

            if self._sf_.TW.func.body.__class__.__name__ == '_2dCSCG_ScalarField':
                #
                if self._sf_.func.ftype == 'standard':
                    return self._standard_scalar_(update_cochain=update_cochain)
                else:
                    raise NotImplementedError(f"2dCSCG outer 0-form cannot (target func) "
                                              f"discretize _2dCSCG_ScalarField of ftype={self._sf_.func.ftype}")

            else:
                raise NotImplementedError(f'2dCSCG outer 0-form can not (target func) '
                                          f'discretize {self._sf_.TW.func.body.__class__}.')

        elif target == 'BC':
                raise NotImplementedError(f'2dCSCG outer 0-form can not (target BC) '
                                          f'discretize {self._sf_.TW.BC.body.__class__}.')


        else:
            raise NotImplementedError(f"2dCSCG outer 0-form cannot discretize "
                                      f"while targeting at {target}.")