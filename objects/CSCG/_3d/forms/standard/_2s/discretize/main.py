# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly

from objects.CSCG._3d.forms.standard._2s.discretize.vector.boundary_wise import _3dCSCG_Discretize_BoundaryWise
from objects.CSCG._3d.forms.standard._2s.discretize.vector.standard import _3dCSCG_Discretize_Standard


class _3dCSCG_Discretize(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._standard_ = _3dCSCG_Discretize_Standard(sf)
        self._boundary_wise_ = _3dCSCG_Discretize_BoundaryWise(sf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', **kwargs):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :param kwargs:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        SELF = self._sf_

        if target == 'func':

            if SELF.CF.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.CF.ftype == 'standard':
                    return self._standard_(update_cochain=update_cochain, **kwargs)

                else:
                    raise NotImplementedError(f"3dCSCG 2-form cannot (target func) "
                                              f"discretize _3dCSCG_VectorField of ftype={SELF.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 2-form can not (target func) '
                                          f'discretize {SELF.CF.__class__}.')

        elif target == 'BC':
            if SELF.BC.CF.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.BC.CF.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self._standard_(update_cochain=False, target='BC', **kwargs)

                elif SELF.BC.CF.ftype == "boundary-wise":
                    # we will always not update cochain & and always set target to be "BC"
                    return self._boundary_wise_(**kwargs)

                else:
                    raise NotImplementedError(f"3dCSCG 2-form cannot (target BC) "
                                              f"discretize _3dCSCG_VectorField of ftype={SELF.BC.CF.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 2-form can not (target BC) '
                                          f'discretize {SELF.BC.CF.__class__}.')

        else:
            raise NotImplementedError(f"3dCSCG 2-form cannot discretize while targeting at {target}.")