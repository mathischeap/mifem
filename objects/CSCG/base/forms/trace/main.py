# -*- coding: utf-8 -*-
from objects.CSCG.base.forms.trace.whether import CSCG_trace_form_Whether
from objects.CSCG.base.forms.trace.num import CSCG_trace_form_NUM


# noinspection PyUnresolvedReferences
class CSCG_Trace_Form:
    """"""
    def ___init___(self):
        self._whether_ = CSCG_trace_form_Whether(self)
        self._num_ = CSCG_trace_form_NUM(self)

    @property
    def orientation(self):
        """(str) The orientation."""
        return self._orientation_

    @property
    def whether(self):
        return self._whether_

    @property
    def num(self):
        return self._num_