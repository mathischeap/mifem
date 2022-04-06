
from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard._2s.main import _3dCSCG_2Form
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from objects.CSCG._3d.forms.standard._1s.special.vortex_detection import ___3dCSCG_1Form_Vortex_Detection___
from objects.CSCG._3d.forms.standard._1s.special.helpers.cross_product_1__ip_1 import ___3dCSCG_1Form_CrossProduct_1__ip_1___
from objects.CSCG._3d.forms.standard._1s.special.helpers.cross_product_2__ip_2 import ___3dCSCG_1Form_CrossProduct_2__ip_2___




class _1Form_Special(FrozenOnly):
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product_1f__ip_1f(self, u, e, quad_degree=None, output='2-M-1'):
        """
        (self X 1-form, 1-form). To first cross product with a 1-form then do a inner product with
        another 1-form.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e(column) and u (row)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_1Form_CrossProduct_1__ip_1___(self._sf_, u, e, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    def cross_product_2f__ip_2f(self, u, e, quad_degree=None, output='2-M-1'):
        """
        (self X 2-form, 2-form). To first cross product with a 2-form then do an inner product with
        another 2-form.

        output:
            '2-M-1': Means we return a local matrix refers to local dofs of e (column) and u (row)

        :return:
        """
        if output == '2-M-1':
            SCP_generator = ___3dCSCG_1Form_CrossProduct_2__ip_2___(self._sf_, u, e, quad_degree=quad_degree)
        else:
            raise NotImplementedError(f"output={output} is not implemented.")

        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')




    def ___PRIVATE_projected_into_2form_exactly___(self):
        """We project this 1form into a 2form exactly. Since it is an
        exact projection, we will use a space one degree higher
        than the space of this 1form. The mesh will be the same mesh.
        """
        space = self._sf_.space
        mesh = self._sf_.mesh

        sp = space.p
        op = list()
        for i in sp: op.append(i+1)

        SPACE = space.__class__(op, None)

        f2 = _3dCSCG_2Form(mesh, SPACE,
                           is_hybrid=self._sf_.IS.hybrid,
                           orientation=self._sf_.orientation,
                           numbering_parameters=self._sf_.numbering._numbering_parameters_,
                           name='Projected_2form_of_'+self._sf_.standard_properties.name
                           )

        W21 = self._sf_.operators.wedge(f2)
        invM2 = f2.matrices.mass.inv

        lc1 = self._sf_.cochain.local

        lc2 = dict()
        for i in lc1:
            lc2[i] = invM2[i] @ W21[i] @ lc1[i]

        f2.cochain.local = lc2

        return f2

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_1Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_