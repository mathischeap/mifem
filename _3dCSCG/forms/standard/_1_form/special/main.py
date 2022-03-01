
from screws.frozen import FrozenOnly
import numpy as np
from _3dCSCG.forms.standard._2_form.main import _3dCSCG_2Form
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from scipy import sparse as spspa
from _3dCSCG.forms.standard._1_form.special.vortex_detection import ___3dCSCG_1Form_Vortex_Detection___



class _1Form_Special(FrozenOnly):
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product(self, u, e, quad_degree=None):
        """
        We do ``(self X other, e)`` where ``self`` and ``other`` both are n-form, n be either 1 or 2.

        :return:
        """
        SCP_generator = ___3dCSCG_1Form_CrossProduct___(self._sf_, u, e, quad_degree=quad_degree)
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

        f2 = _3dCSCG_2Form(mesh, SPACE, is_hybrid=self._sf_.IS_hybrid,
                           orientation=self._sf_.orientation,
                           numbering_parameters=self._sf_.numbering._numbering_parameters_,
                           name='Projected_2form_of_'+self._sf_.standard_properties.name)


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


class ___3dCSCG_1Form_CrossProduct___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(\\omega \\times u, e)`.

    :param w1:
    :param u1:
    :param e1:
    :param quad_degree:
    """
    def __init__(self, w1, u1, e1, quad_degree=None):
        assert u1.ndim == w1.ndim == e1.ndim, " <_3dCSCG_1Form_CrossProduct> "
        assert u1.k == w1.k == e1.k == 1, " <_3dCSCG_1Form_CrossProduct> "
        assert u1.mesh == w1.mesh, "Meshes do not match."
        assert u1.mesh == e1.mesh, "Meshes do not match."
        self._mesh_ = u1.mesh
        self._u_ = u1
        self._w_ = w1
        self._e_ = e1

        if quad_degree is None:
            quad_degree = [int(np.max([u1.dqp[i], w1.dqp[i], e1.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = u1.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, wbf = w1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=True)
        _, ubf = u1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        if e1 is u1:
            ebf = ubf
        else:
            _, ebf = e1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._xietasigma_ = xietasigma
        self._qw_ = quad_weights
        self._wbf_ = wbf
        self._ubf_ = ubf
        self._ebf_ = ebf
        self._JM_ = self._mesh_.elements.coordinate_transformation.Jacobian_matrix(*xietasigma)
        self._iJ_ = self._mesh_.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=self._JM_)
        self.DO_reset_cache()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._J_cache_ = dict()

    def _J_(self, i):
        element = self._mesh_.elements[i]
        typeWr2Metric = element.type_wrt_metric.mark
        if typeWr2Metric in self._J_cache_:
            return self._J_cache_[typeWr2Metric]
        else:
            JM = self._JM_[i]
            if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
                J00 = JM[1][1] * JM[2][2]
                J01 = 0
                J02 = 0
                J10 = 0
                J11 = JM[2][2] * JM[0][0]
                J12 = 0
                J20 = 0
                J21 = 0
                J22 = JM[0][0] * JM[1][1]
            else:
                J00 = JM[1][1] * JM[2][2] - JM[1][2] * JM[2][1]
                J01 = JM[2][1] * JM[0][2] - JM[2][2] * JM[0][1]
                J02 = JM[0][1] * JM[1][2] - JM[0][2] * JM[1][1]
                J10 = JM[1][2] * JM[2][0] - JM[1][0] * JM[2][2]
                J11 = JM[2][2] * JM[0][0] - JM[2][0] * JM[0][2]
                J12 = JM[0][2] * JM[1][0] - JM[0][0] * JM[1][2]
                J20 = JM[1][0] * JM[2][1] - JM[1][1] * JM[2][0]
                J21 = JM[2][0] * JM[0][1] - JM[2][1] * JM[0][0]
                J22 = JM[0][0] * JM[1][1] - JM[0][1] * JM[1][0]
            J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
            # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
            self._J_cache_[typeWr2Metric] = J
            return J

    def __call__(self, i):
        typeWr2Metric = self._mesh_.elements[i].type_wrt_metric.mark

        u0, u1, u2 =  self._ubf_ # a
        w0, w1, w2 =  self._wbf_ # b; given
        e0, e1, e2 =  self._ebf_ # epsilon

        b0p = np.einsum('ij, i -> j', w0, self._w_.cochain.local_('x')[i], optimize='greedy')
        b1p = np.einsum('ij, i -> j', w1, self._w_.cochain.local_('y')[i], optimize='greedy')
        b2p = np.einsum('ij, i -> j', w2, self._w_.cochain.local_('z')[i], optimize='greedy')

        iJ = self._iJ_[i]
        J00, J01, J02, J10, J11, J12, J20, J21, J22 = self._J_(i)

        if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
            b0 = b0p * iJ[0][0]
            b1 = b1p * iJ[1][1]
            b2 = b2p * iJ[2][2]

            B01, B02 = -iJ[0][0]*b2, iJ[0][0]*b1
            B10, B12 = iJ[1][1]*b2, -iJ[1][1]*b0
            B20, B21 = -iJ[2][2]*b1, iJ[2][2]*b0

            m01 = B10*J00
            m02 = B20*J00
            m10 = B01*J11
            m12 = B21*J11
            m20 = B02*J22
            m21 = B12*J22

            # put `-` because a x b = - b x a
            M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
            M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')

            M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
            M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')

            M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
            M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')

            M = ([None, spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
                 [spspa.csc_matrix(M10), None, spspa.csc_matrix(M12)],
                 [spspa.csc_matrix(M20), spspa.csc_matrix(M21), None])
        else:
            b0 = b0p * iJ[0][0] + b1p * iJ[1][0] + b2p * iJ[2][0]
            b1 = b0p * iJ[0][1] + b1p * iJ[1][1] + b2p * iJ[2][1]
            b2 = b0p * iJ[0][2] + b1p * iJ[1][2] + b2p * iJ[2][2]

            B00, B01, B02 = (iJ[0][1]*b2-iJ[0][2]*b1), (iJ[0][2]*b0-iJ[0][0]*b2), (iJ[0][0]*b1-iJ[0][1]*b0)
            B10, B11, B12 = (iJ[1][1]*b2-iJ[1][2]*b1), (iJ[1][2]*b0-iJ[1][0]*b2), (iJ[1][0]*b1-iJ[1][1]*b0)
            B20, B21, B22 = (iJ[2][1]*b2-iJ[2][2]*b1), (iJ[2][2]*b0-iJ[2][0]*b2), (iJ[2][0]*b1-iJ[2][1]*b0)

            m00 = B00*J00 + B01*J01 + B02*J02
            m01 = B10*J00 + B11*J01 + B12*J02
            m02 = B20*J00 + B21*J01 + B22*J02
            m10 = B00*J10 + B01*J11 + B02*J12
            m11 = B10*J10 + B11*J11 + B12*J12
            m12 = B20*J10 + B21*J11 + B22*J12
            m20 = B00*J20 + B01*J21 + B02*J22
            m21 = B10*J20 + B11*J21 + B12*J22
            m22 = B20*J20 + B21*J21 + B22*J22

            # put `-` because a x b = - b x a
            M00 = - np.einsum('iw, jw, w -> ij', e0, u0, m00*self._qw_, optimize='greedy')
            M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
            M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')

            M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
            M11 = - np.einsum('iw, jw, w -> ij', e1, u1, m11*self._qw_, optimize='greedy')
            M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')

            M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
            M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')
            M22 = - np.einsum('iw, jw, w -> ij', e2, u2, m22*self._qw_, optimize='greedy')

            M = ([spspa.csc_matrix(M00), spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
                 [spspa.csc_matrix(M10), spspa.csc_matrix(M11), spspa.csc_matrix(M12)],
                 [spspa.csc_matrix(M20), spspa.csc_matrix(M21), spspa.csc_matrix(M22)])
        MW = spspa.bmat(M, format='csc')

        return MW

