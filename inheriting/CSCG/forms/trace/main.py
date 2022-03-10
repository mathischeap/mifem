
from inheriting.CSCG.forms.trace.IS import CSCG_trace_form_IS
from inheriting.CSCG.forms.trace.num import CSCG_trace_form_NUM


# noinspection PyUnresolvedReferences
class CSCG_Trace_Form:
    """"""
    def ___init___(self):
        self._IS_ = CSCG_trace_form_IS(self)
        self._num_ = CSCG_trace_form_NUM(self)

    @property
    def GLOBAL_num_dofs(self):
        """(int) Return the total number of dofs this trace form has."""
        return self.numbering.gathering.GLOBAL_num_dofs

    @property
    def orientation(self):
        """(str) The orientation."""
        return self._orientation_

    @property
    def IS(self):
        return self._IS_

    @property
    def num(self):
        return self._num_