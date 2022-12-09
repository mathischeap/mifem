# -*- coding: utf-8 -*-
""""""
from objects.CSCG.base.forms.standard.whether import CSCG_standard_form_Whether
from objects.CSCG.base.forms.standard.num import CSCG_standard_form_NUM

# noinspection PyUnresolvedReferences
class CSCG_Standard_Form:
    """"""
    def ___init___(self):
        self._whether_ = CSCG_standard_form_Whether(self)
        self._num_ = CSCG_standard_form_NUM(self)

    @property
    def orientation(self):
        return self._orientation_

    @property
    def whether(self):
        return self._whether_

    @property
    def num(self):
        return self._num_