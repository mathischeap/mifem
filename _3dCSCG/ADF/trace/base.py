""" """

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.frozen import FrozenOnly
from _3dCSCG.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE

from INHERITING.CSCG.ADF.trace.main_BASE import CSCG_Algebra_DUAL_Trace_Form
from scipy import sparse as spspa

from TOOLS.linear_algebra.elementwise_cache import EWC_SparseMatrix


class _3dCSCG_Algebra_DUAL_Trace_Form(CSCG_Algebra_DUAL_Trace_Form, _3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""
    def __init__(self, ndim, mesh, space, orientation, name):
        """"""
        super().__init__(ndim, mesh, space)
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_trace')
        self._DO_ = _3dCSCG_Algebra_DUAL_Trace_Form_DO(self)
        self._cochain_ = _3dCSCG_Algebra_DUAL_Trace_Form_Cochain(self)
        self._coboundary_ = _3dCSCG_Algebra_DUAL_Trace_Form_Coboundary(self)

    def RESET_cache(self):
        """"""

    def ___PreFrozenChecker___(self):
        """This method will be run automatically before we freeze the object. This is very important because
        we have to make sure that the information is consistent across the prime and the dual.
        """
        assert self.prime.k == self.k           , "Trivial check k."
        assert self.prime.mesh == self.mesh     , "Trivial check mesh."
        assert self.prime.space == self.space   , "Trivial check space."
        assert self.prime.orientation == self.orientation, "Trivial check orientation."
        assert self.prime.IS_hybrid is self.IS_hybrid    , "prime must be hybrid form, now it is not."

    def ___PRIVATE_generate_mass_matrix___(self):
        """For algebra dual forms, this method will only be called once. The result (a mass matrix) will be cached in
        `self._mass_matrix_` because we think it is an essential property for algebra dual forms and it will be used
        for multiple times. Therefore, we cache it somewhere.

        :return: A tuple of two outputs: the mass matrix and the inverse mass matrix.
        """
        TEW_mass = self.prime.matrices.mass
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
            self.mesh.elements.___PRIVATE_elementwise_cache_metric_key___, # cache keys as mesh element types wrt metric.
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




class _3dCSCG_Algebra_DUAL_Trace_Form_DO(FrozenOnly):
    def __init__(self, dt):
        self._dt_ = dt
        self._freeze_self_()




class _3dCSCG_Algebra_DUAL_Trace_Form_Cochain(FrozenOnly):
    """The cochain of algebra dual form is equal to the mass matrix dot the cochain of the prime form.

    dual_cochain = mass_matrix dot prime_cochain.

    """
    def __init__(self, dt):
        """
        :param dt: The dual standard form.
        """
        self._dt_ = dt
        self._local_ = None  # this is a key property, should not reset it.
        self._freeze_self_()

    @property
    def local(self):
        """We know that the local cochain of the prime form is a dict whose keys are local element numbers and values
        are the local cochains (1-d array). While for algebra dual standard forms, we make a EWC_ColumnVector for it
        since we do not want to save the local cochain of the algebra dual form. We will generate the cochain when
        we call it in real time.

        :return:
        """
        if self._local_ is None:
            self._local_ = ____3dCSCG_ADTF_Cochain_Local____(self)
            # the local cochain will be renew automatically if the local cochain of the prime form is renewed.
        return self._local_

    @local.setter
    def local(self, local):
        """"""
        assert len(local) == len(self._dt_.mesh.elements), \
            f"length of local is wrong!"
        prime_local_cochain = dict()
        iM = self._dt_.inverse_mass_matrix
        for i in self._dt_.mesh.elements:
            assert i in local, \
                f"mesh element #{i} is missing in local in Core #{rAnk}."
            prime_local_cochain[i] = iM[i] @ local[i]
        self._dt_.prime.cochain.local = prime_local_cochain

    @property
    def globe(self):
        raise NotImplementedError()

    @globe.setter
    def globe(self, globe):
        """We have to set the cochain to the local cochain of the prime.

        :param globe:
        :return:
        """
        raise NotImplementedError()

    def __getitem__(self, i):
        """If `i` is an element number in this core, we should be able to return a local cochain of it (if not None)"""
        return self.local[i]

    def __contains__(self, i):
        """Return if an element number (`i`) is in this core."""
        return i in self.local

    def __iter__(self):
        """Go through all mesh element numbers in this core."""
        for i in self.local:
            yield i

    def __len__(self):
        """Actually return how many mesh elements in this core."""
        return len(self.local)




class ____3dCSCG_ADTF_Cochain_Local____(FrozenOnly):
    """"""
    def __init__(self, dt_CO):
        """"""
        self._PC_ = dt_CO._dt_.prime.cochain
        self._MM_ = dt_CO._dt_.mass_matrix
        self._freeze_self_()

    def __getitem__(self, i):
        """"""
        return self._MM_[i] @ self._PC_.local[i]

    def __contains__(self, i):
        """"""
        return i in self._PC_.local

    def __iter__(self):
        """Go through all mesh element numbers in this core."""
        for i in self._PC_.local:
            yield i

    def __len__(self):
        """Actually return how many mesh elements in this core."""
        return len(self._PC_.local)




class _3dCSCG_Algebra_DUAL_Trace_Form_Coboundary(FrozenOnly):
    """This is one of the key properties of a standard algebraic dual form. To perform the coboundary, now we need the
    help of another dual trace form.
    """
    def __init__(self, dt):
        self._dt_ = dt
        self._T_ = None
        self._freeze_self_()

    @property
    def trace_matrix(self):
        """"""
        if self._T_ is None:
            if self._dt_.k == 2:
                formName = '_2Trace'
            elif self._dt_.k == 1:
                formName = '_1Trace'
            elif self._dt_.k == 0:
                formName = '_0Trace'
            else:
                raise Exception()
            T = getattr(self._dt_.space.trace_matrix, formName)[0].T
            self._T_ = \
                EWC_SparseMatrix(self._dt_.mesh.elements, T, 'constant')
        return self._T_








if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\trace\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))([12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)