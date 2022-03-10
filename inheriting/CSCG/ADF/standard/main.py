


from inheriting.CSCG.ADF.standard.IS import CSCG_ADF_SF_IS


# noinspection PyUnresolvedReferences
class CSCG_Algebra_DUAL_Standard_Form:
    """"""
    def ___init___(self):
        self._IS_ = CSCG_ADF_SF_IS(self)

    @property
    def orientation(self):
        """An AD standard form can be either inner or outer."""
        return self._orientation_


    @property
    def IS(self):
        return self._IS_