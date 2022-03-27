

from screws.freeze.base import FrozenOnly
from inheriting.CSCG.forms.base.time_wise.func.do import CSCG_FTWF_DO


class CSCG_Form_TimeWise_Func(FrozenOnly):
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWF_DO(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        """Clear cache."""
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None

    def ___DO_set_func_body_as___(self, scalar_vector_or_ES, variable_name=None):
        """
        We use this function such that we can set func from exact_solution object.

        If we want the func can be restored, we have to use this function, and we only store and restore
        func from exact_solution object.

        :param scalar_vector_or_ES: scalar, vector or exact_solution.
        :param variable_name:
        :return:
        """
        if scalar_vector_or_ES.__class__.__name__ in (f'_{self._ndim_}dCSCG_ScalarField',
                                                      f'_{self._ndim_}dCSCG_VectorField',
                                                      f'_{self._ndim_}dCSCG_TensorField'):
            assert variable_name is None
            self.body = scalar_vector_or_ES  # has to pass the checker.
        else:
            # only by setting ES and variable_name, the func can be restored.
            self_body = getattr(scalar_vector_or_ES.status, variable_name)
            self._f_.___PRIVATE_TW_FUNC_body_checker___(self_body)  # has to pass the checker.
            # do not use body.setter because we want to set _ES_ and _ES_variable_name_ attributes.
            self._body_ = self_body
            self._ES_ = scalar_vector_or_ES
            self._ES_variable_name_ = variable_name

    @property
    def body(self):
        """The function body."""
        return self._body_

    @body.setter
    def body(self, body):
        self._f_.___PRIVATE_TW_FUNC_body_checker___(body)  # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The type of the function body: `self.body`."""
        return self.body.ftype

    @property
    def do(self):
        return self._DO_

