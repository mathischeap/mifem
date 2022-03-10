

from inheriting.CSCG.ADF.trace.IS import CSCG_ADT_TF_IS


# noinspection PyUnresolvedReferences
class CSCG_Algebra_DUAL_Trace_Form:
    """"""
    def ___init___(self):
        self._IS_ = CSCG_ADT_TF_IS(self)

    @property
    def orientation(self):
        return self._orientation_


    @property
    def IS(self):
        return self._IS_