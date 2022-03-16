



from _3dCSCG.forms.standard._2_form.special.vortex_detection import ___3dCSCG_2Form_Vortex_Detection___

from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix

from _3dCSCG.forms.standard._2_form.special.helpers.cross_product_2__ip_2 import ___3dCSCG_2Form_CrossProduct_2__ip_2___




class _2Form_Special(FrozenOnly):
    def __init__(self, _2sf):
        self._sf_ = _2sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product_2f__ip_2f(self, u, e, quad_degree=None):
        """
        (self X 2form, 2form)

        (self X u    , e    )

        We do ``(self X other, e)`` where ``self`` and ``other`` both are n-form, n be either 1 or 2.

        :return:
        """
        SCP_generator = ___3dCSCG_2Form_CrossProduct_2__ip_2___(self._sf_, u, e, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_2Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_



