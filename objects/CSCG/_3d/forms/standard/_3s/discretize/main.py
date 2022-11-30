from components.freeze.base import FrozenOnly

from objects.CSCG._3d.forms.standard._3s.discretize.scalar.standard import _3dCSCG_Discretize_Standard


class _3dCSCG_Discretize(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._standard_ = _3dCSCG_Discretize_Standard(sf)
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
            if SELF.CF.__class__.__name__ == '_3dCSCG_ScalarField':

                if SELF.CF.ftype == 'standard':
                    return self._standard_(update_cochain=update_cochain, **kwargs)
                else:
                    raise NotImplementedError(
                        f"3dCSCG 3-form cannot (target func) discretize _3dCSCG_ScalarField of ftype={SELF.CF.ftype}")

            else:
                raise NotImplementedError(
                    f'3dCSCG 3-form can not (target func) discretize {SELF.CF.__class__}.')
        else:
            raise NotImplementedError(f"3dCSCG 3-form cannot discretize while targeting at {target}.")