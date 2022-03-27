from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard._1_form.outer.discretize.vector.standard import _2dCSCG_S1Fo_Discretize_StandardVector

class _2dCSCG_S1Fo_Discretize(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._standard_vector_ = _2dCSCG_S1Fo_Discretize_StandardVector(sf)
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', **kwargs):
        """
        Discretize the current function (a vector field:
        :class:`_3dCSCG.form.continuous.vector._3dCSCG_VectorField`) to cochain.
        It is actually a wrapper of multiple methods that discretize functions of different types (a vector
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':


            if self._sf_.TW.func.body.__class__.__name__ == '_2dCSCG_VectorField':

                if self._sf_.func.ftype == 'standard':
                    return self._standard_vector_(
                        update_cochain=update_cochain, target='func', **kwargs)
                else:
                    raise NotImplementedError()

            else:
                raise Exception()


        elif target == 'BC':
            raise NotImplementedError(f'2dCSCG 1-form can not (target BC) '
                                      f'discretize {self._sf_.TW.BC.body.__class__}.')


        else:
            raise NotImplementedError(f"2dCSCG 1-form cannot discretize "
                                      f"while targeting at {target}.")