
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.Tr._2Tr.discretize.vector.standard import _3dCSCG_2Tr_Discretize_StandardVector
from objects.CSCG._3d.forms.Tr._2Tr.discretize.vector.boundary_wise import _3dCSCG_2Tr_Discretize_BoundaryWiseVector

from objects.CSCG._3d.forms.Tr._2Tr.discretize.scalar.standard import _3dCSCG_2Tr_Discretize_StandardScalar
from objects.CSCG._3d.forms.Tr._2Tr.discretize.scalar.boundary_wise import _3dCSCG_2Tr_Discretize_BoundaryWiseScalar




class _3dCSCG_2Tr_Discretize(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._standard_vector_ = _3dCSCG_2Tr_Discretize_StandardVector(Tr)
        self._boundary_wise_vector_ = _3dCSCG_2Tr_Discretize_BoundaryWiseVector(Tr)
        self._standard_scalar_ = _3dCSCG_2Tr_Discretize_StandardScalar(Tr)
        self._boundary_wise_scalar_ = _3dCSCG_2Tr_Discretize_BoundaryWiseScalar(Tr)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param target:
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """
        SELF = self._Tr_

        if target == 'func':
            if SELF.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.func.ftype == 'standard':
                    return self._standard_scalar_(
                        update_cochain=update_cochain, **kwargs)
                else:
                    raise Exception(f'3dCSCG 2-Tr can not (target func) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.func.ftype}.')

            elif SELF.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':
                if SELF.func.ftype == 'standard': # we will discretize the norm component of the vector.
                    return self._standard_vector_(
                        update_cochain=update_cochain, **kwargs)
                else:
                    raise Exception(f'3dCSCG 2-Tr can not (target func) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.func.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 2-Tr can not (target func) '
                                          f'discretize {SELF.TW.func.body.__class__}.')

        elif target == 'BC': # We target at the BC, so we do not update the cochain!

            if SELF.TW.BC.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if SELF.BC.ftype == 'standard':
                    return self._standard_scalar_(
                        update_cochain=False, target='BC', **kwargs)
                elif SELF.BC.ftype == 'boundary-wise':
                    return self._boundary_wise_scalar_(
                        **kwargs) # must be False update_cochain and 'BC' target.
                else:
                    raise Exception(f'3dCSCG 2-Tr can not (target BC) discretize '
                                    f'_3dCSCG_ScalarField of ftype {SELF.BC.ftype}.')

            elif SELF.TW.BC.body.__class__.__name__ == '_3dCSCG_VectorField':
                if SELF.BC.ftype == 'standard': # we will discretize the norm flux of the vector.
                    return self._standard_vector_(
                        update_cochain=False, target='BC', **kwargs)
                elif SELF.BC.ftype == 'boundary-wise': # we will discretize the norm flux of the vector.
                    return self._boundary_wise_vector_(
                        **kwargs) # must be False update_cochain and 'BC' target.
                else:
                    raise Exception(f'3dCSCG 2-Tr can not (target BC) discretize '
                                    f'_3dCSCG_VectorField of ftype {SELF.BC.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 2-Tr can not (target BC) '
                                          f'discretize {SELF.TW.BC.body.__class__}.')

        else:
            raise NotImplementedError(f"target={target} not implemented "
                                      f"for 3d CSCG 2-Tr form discretization.")



if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_2_trace\discretize\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)