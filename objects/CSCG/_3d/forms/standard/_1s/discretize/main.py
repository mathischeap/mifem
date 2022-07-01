# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly

from objects.CSCG._3d.forms.standard._1s.discretize.vector.boundary_wise import _3dCSCG_Discretize_BoundaryWise
from objects.CSCG._3d.forms.standard._1s.discretize.vector.standard import _3dCSCG_Discretize_Standard


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
            if SELF.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.func.ftype == 'standard':
                    return self._standard_(update_cochain=update_cochain, **kwargs)
                else:
                    raise NotImplementedError(f"3dCSCG 1-form cannot (target func) discretize "
                                              f"_3dCSCG_VectorField of ftype={SELF.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 1-form can not discretize {SELF.TW.func.body.__class__}.')

        elif target == 'BC':
            if SELF.TW.BC.body.__class__.__name__ == '_3dCSCG_VectorField':

                if SELF.BC.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self._standard_(update_cochain=False, target='BC', **kwargs)
        
        else:
            raise NotImplementedError(f"3dCSCG 1-form cannot discretize while targeting at {target}.")