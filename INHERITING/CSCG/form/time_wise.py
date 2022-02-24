# -*- coding: utf-8 -*-

"""
If one type of form (like standard form, trace form or whatever) want to have its ``func`` mode, it has
to inherit ``CSCG_Form_Func`` and ``CSCG_Form_TimeWise``.
"""

from SCREWS.frozen import FrozenOnly



class CSCG_Form_TimeWise(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self._current_time_ = None

        self._func_ = CSCG_Form_TimeWise_Func(self)
        self._BC_ = CSCG_Form_TimeWise_BC(self)

        self._DO_ = CSCG_Form_TimeWise_DO(self)
        self._freeze_self_()

    @property
    def current_time(self):
        assert self._current_time_ is not None, "current_time is None, set it firstly."
        return self._current_time_

    @current_time.setter
    def current_time(self, current_time):
        assert isinstance(current_time, (int, float))
        self._current_time_ = current_time




    @property
    def func(self):
        return self._func_

    @property
    def BC(self):
        """Boundary condition."""
        return self._BC_




    @property
    def DO(self):
        return self._DO_


    def ___DO_push_func_to_instant___(self, t=None):
        """Will update current_time if t is not None"""
        if t is not None: self.current_time = t
        self._f_.func._body_ = self.func.body.___DO_evaluate_func_at_time___(self.current_time)
        self._f_.func._ftype_ = self.func.body.ftype


    def ___DO_push_BC_to_instant___(self, t=None):
        """Will update current_time if t is not None"""
        if t is not None: self.current_time = t
        self._f_.BC._partial_cochain_ = None # IMPORTANT: we clean the previous PartialCochain instance.
        self._f_.BC._body_ = self.BC.body.___DO_evaluate_func_at_time___(self.current_time)
        self._f_.BC._ftype_ = self.BC.body.ftype


    def ___DO_push_all_to_instant___(self, t=None):
        """Will update current_time if t is not None"""
        if self.func.body is not None:
            self.___DO_push_func_to_instant___(t=t)
        if self.BC.body is not None:
            self.___DO_push_BC_to_instant___(t=t)




class CSCG_Form_TimeWise_DO(FrozenOnly):
    def __init__(self, tw):
        self._tw_ = tw
        self._freeze_self_()

    def push_func_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_func_to_instant___(t=t)

    def push_BC_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_BC_to_instant___(t=t)

    def push_all_to_instant(self, t=None):
        """Will update current_time if t is not None"""
        self._tw_.___DO_push_all_to_instant___(t=t)






# ---------- FUNC ------------------------------------------------------------------------------------------------------

class CSCG_Form_TimeWise_Func(FrozenOnly):
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWF_DO(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
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
            self._f_.___TW_FUNC_body_checker___(self_body)  # has to pass the checker.
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
        self._f_.___TW_FUNC_body_checker___(body)  # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The type of the function body: `self.body`."""
        return self.body.ftype

    @property
    def DO(self):
        return self._DO_


class CSCG_FTWF_DO(FrozenOnly):
    def __init__(self, twf):
        self._twf_ = twf
        self._freeze_self_()

    def set_func_body_as(self, *args, **kwargs):
        self._twf_.___DO_set_func_body_as___(*args, **kwargs)









class CSCG_Form_Func(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self._body_ = None
        self._ftype_ = None

    @property
    def body(self):
        """
        The function body.
        """
        return self._body_

    @property
    def ftype(self):
        """
        The function type.

        :return: The function type.
        :rtype: str
        """
        return self._ftype_






# ------------------- BC --------------------------------------------------------------------------------

class CSCG_Form_TimeWise_BC(FrozenOnly):
    """Note that the TW.BC of a form will not be saved!!!!.
    """
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWBC_DO(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        """Clear cache."""
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None

    def ___DO_set_BC_body_as___(self, scalar_vector):
        """
        We use this function such that we can set BC.

        :param scalar_vector: scalar, vector
        :return:
        """
        self.body = scalar_vector  # has to pass the checker.


    @property
    def body(self):
        """The function body."""
        return self._body_

    @body.setter
    def body(self, body):
        self._f_.___TW_BC_body_checker___(body)  # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The function type."""
        return self.body.ftype

    @property
    def DO(self):
        return self._DO_


class CSCG_FTWBC_DO(FrozenOnly):
    def __init__(self, twBC):
        self._twBC_ = twBC
        self._freeze_self_()

    def set_BC_body_as(self, *args, **kwargs):
        self._twBC_.___DO_set_BC_body_as___(*args, **kwargs)









class CSCG_Form_BC_func(FrozenOnly):
    def __init__(self, f):
        self._f_ = f

    def RESET_cache(self):
        self._body_ = None
        self._ftype_ = None

    @property
    def body(self):
        """
        The function body.
        """
        return self._body_

    @property
    def ftype(self):
        """
        The function type.

        :return: The function type.
        :rtype: str
        """
        return self._ftype_





#===================================== ABOVE ===================================================================