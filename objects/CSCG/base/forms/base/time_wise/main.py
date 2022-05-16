# -*- coding: utf-8 -*-

"""
If one type of form (like standard form, trace form or whatever) want to have its ``func`` mode, it has
to inherit ``CSCG_Form_Func`` and ``CSCG_Form_TimeWise``.
"""

from screws.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.time_wise.BC.main import CSCG_Form_TimeWise_BC
from objects.CSCG.base.forms.base.time_wise.func.main import CSCG_Form_TimeWise_Func
from objects.CSCG.base.forms.base.time_wise.do import CSCG_Form_TimeWise_DO


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
    def do(self):
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
