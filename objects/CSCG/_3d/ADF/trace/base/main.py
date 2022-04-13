""" """

import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE

from scipy import sparse as spspa

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix

from objects.CSCG._3d.ADF.trace.base.do import _3dCSCG_Algebra_DUAL_Trace_Form_DO
from objects.CSCG._3d.ADF.trace.base.cochain.main import _3dCSCG_Algebra_DUAL_Trace_Form_Cochain
from objects.CSCG._3d.ADF.trace.base.coboundary import _3dCSCG_Algebra_DUAL_Trace_Form_Coboundary
from objects.CSCG._3d.ADF.trace.base.IS import _3dCSCG_ADT_TF_IS
from objects.CSCG._3d.ADF.trace.base.num import _3dCSCG_ADF_T_NUM
from objects.CSCG._3d.ADF.trace.base.matrices import _3dCSCG_ADT_TF_Matrices




class _3dCSCG_Algebra_DUAL_Trace_Form(_3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""
    def __init__(self, ndim, mesh, space, prime, orientation, name):
        """"""
        super(_3dCSCG_Algebra_DUAL_Trace_Form, self).__init__(ndim, mesh, space, prime)

        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_trace')
        self._DO_ = _3dCSCG_Algebra_DUAL_Trace_Form_DO(self)
        self._cochain_ = _3dCSCG_Algebra_DUAL_Trace_Form_Cochain(self)
        self._coboundary_ = _3dCSCG_Algebra_DUAL_Trace_Form_Coboundary(self)
        self._IS_ = None
        self._num_ = None
        self._matrices_ = None

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

    def ___PRIVATE_generate_mass_matrix___(self):
        """For algebra dual forms, this method will only be called once.
        The result (a mass matrix) will be cached in
        `self._mass_matrix_` because we think it is an essential property for algebra dual forms,
        and it will be used
        for multiple times. Therefore, we cache it somewhere.

        :return: A tuple of two outputs: the mass matrix and the inverse mass matrix.
        """
        TEW_mass = self.prime.matrices.mass_TEW
        # we now need to make a mesh-element-wise mass matrix
        MAP = self.mesh.trace.elements.map
        local_cache = dict()
        MEW_mass = dict()

        for i in MAP:
            element = self.mesh.elements[i]
            mark = element.type_wrt_metric.mark

            if isinstance(mark, str) and mark in local_cache:
                MEW_mass[i] = local_cache[mark]
            else:

                M = list()
                for t in MAP[i]:
                    assert t in TEW_mass, "A trivial check!"
                    Mt = TEW_mass[t]
                    M.append(Mt)

                M = spspa.bmat([(M[0], None, None, None, None, None),
                                (None, M[1], None, None, None, None),
                                (None, None, M[2], None, None, None),
                                (None, None, None, M[3], None, None),
                                (None, None, None, None, M[4], None),
                                (None, None, None, None, None, M[5])], format='csc')

                if isinstance(mark, str): local_cache[mark] = M

                MEW_mass[i] = M

        M = EWC_SparseMatrix(
            self.mesh.elements, # iterative with mesh elements
            MEW_mass          , # the data as dict
            self.mesh.elements.___PRIVATE_elementwise_cache_metric_key___,
            # cache keys as mesh element types wrt metric.
            )
        iM = M.inv
        return M, iM


    @property
    def DO(self):
        return self._DO_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_ADT_TF_IS(self)
        return self._IS_

    @property
    def num(self):
        if self._num_ is None:
            self._num_ = _3dCSCG_ADF_T_NUM(self)
        return self._num_

    @property
    def matrices(self):
        if self._matrices_ is None:
            self._matrices_ = _3dCSCG_ADT_TF_Matrices(self)
        return self._matrices_





if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\trace\base\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))([12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)