


from screws.freeze.base import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from _2dCSCG.forms.standard._0_form.outer.special.helpers.cross_product_1_i_1 import \
    ___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___



class _0Form_Outer_Special(FrozenOnly):
    def __init__(self, _0sf):
        self._sf_ = _0sf
        self._freeze_self_()


    def cross_product_1f__ip_1f(self, u1, e1, quad_degree=None, output='2-M-1'):
        """
        We do ``(w0 X u1, e1)`` where w0 is self._sf_, u1, e1 are all outer-oriented 1-forms.


        :param u1:
        :param e1:
        :param quad_degree:
        :param output: 0: self, 1: u1, 2: e1.
            '2-M-1': return a mesh-element-wise matrix, row -> e1, col -> u1. So, the cochain
                of self must be known. So it's like we consider `e1` as test functions.
        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___(
                self._sf_, u1, e1, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"not implemented for output={output}.")
        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')