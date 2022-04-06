import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE
from objects.CSCG._3d.ADF.Tr.base.IS import _3dCSCG_ADFTr_form_IS



class _3dCSCG_ADF_Tr_BASE(_3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""
    def __init__(self, ndim, mesh, space, orientation, name):
        super(_3dCSCG_ADF_Tr_BASE, self).__init__(ndim, mesh, space)

        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_Tr')
        self._IS_ = None


    def ___PRIVATE_reset_cache___(self):
        """"""

    def ___PreFrozenChecker___(self):
        """This method will be run automatically before we freeze the object. This is very important because
        we have to make sure that the information is consistent across the prime and the dual.
        """
        assert self.prime.k == self.k           , "Trivial check k."
        assert self.prime.mesh == self.mesh     , "Trivial check mesh."
        assert self.prime.space == self.space   , "Trivial check space."
        assert self.prime.orientation == self.orientation, "Trivial check orientation."
        assert self.prime.IS.hybrid is self.IS.hybrid    , "prime must be hybrid form, now it is not."

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_ADFTr_form_IS(self)
        return self._IS_