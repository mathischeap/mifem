# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.CSCG._2d.forms.standard._0_form.inner.special.helpers.cross_product_1_i_1 import \
    ___2dCSCG_0_i_Form_CrossProduct_0_X_1__ip_1___


class _0Form_Inner_Special(FrozenOnly):
    def __init__(self, _0sf):
        self._sf_ = _0sf
        self._freeze_self_()

    def cross_product_1f__ip_1f(self, u1, e1, quad_degree=None, output='2-M-1'):
        """
        We do ``(w0 X u1, e1)`` where w0 is self._sf_, u1, e1 are all inner oriented 1-forms.

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___2dCSCG_0_i_Form_CrossProduct_0_X_1__ip_1___(
                self._sf_, u1, e1, quad_degree=quad_degree)
        else:
            raise NotImplementedError()

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')