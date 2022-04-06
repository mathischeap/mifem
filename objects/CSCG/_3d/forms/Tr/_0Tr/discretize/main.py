

from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.Tr._0Tr.discretize.vector.standard.flux import _3dCSCG_0Tr_Discretize_StandardVector_Flux
from objects.CSCG._3d.forms.Tr._0Tr.discretize.scalar.standard import _3dCSCG_0Tr_Discretize_StandardScalar


class _3dCSCG_0Tr_Discretize(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._standard_vector_flux_ = _3dCSCG_0Tr_Discretize_StandardVector_Flux(Tr)
        self._standard_scalar_ = _3dCSCG_0Tr_Discretize_StandardScalar(Tr)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param str target:
        :return: The cochain corresponding to the particular discretization scheme.
        """
        SELF = self._Tr_
        if target == 'func':
            if SELF.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.func.ftype == 'standard':
                    return self._standard_scalar_(
                        target = 'func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.func.ftype}.')

            elif SELF.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':
                if SELF.func.ftype == 'standard':
                    return self._standard_vector_flux_(
                        target='func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.func.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target func) '
                                          f'discretize {SELF.TW.func.body.__class__}.')


        elif target == 'BC':
            if SELF.TW.BC.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.BC.ftype == 'standard':
                    return self._standard_scalar_(
                        target = 'BC',
                        update_cochain=False)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target BC) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.BC.ftype}.')
            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target BC) '
                                          f'discretize {SELF.TW.BC.body.__class__}.')
        else:
            raise NotImplementedError()