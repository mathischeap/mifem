
"""

"""
import sys
if './' not in sys.path: sys.path.append('./')


from importlib import import_module
from root.config import *
from SCREWS.frozen import FrozenOnly
from _3dCSCG.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE
from TOOLS.linear_algebra.elementwise_cache import EWC_ColumnVector

from INHERITING.CSCG.ADF.standard.main_BASE import CSCG_Algebra_DUAL_Standard_Form
from scipy.sparse import csr_matrix

from TOOLS.linear_algebra.elementwise_cache import EWC_SparseMatrix





class _3dCSCG_Algebra_DUAL_Standard_Form(CSCG_Algebra_DUAL_Standard_Form, _3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""

    def __init__(self, ndim, mesh, space, orientation, name):
        """"""
        super().__init__(ndim, mesh, space)
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_form')
        self._DO_ = _3dCSCG_Algebra_DUAL_Standard_Form_DO(self)
        self._cochain_ = _3dCSCG_Algebra_DUAL_Standard_Form_Cochain(self)
        self._coboundary_ = _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary(self)

    def RESET_cache(self):
        """"""

    def ___PreFrozenChecker___(self):
        """This method will be run automatically before we freeze the object. This is very important because
        we have to make sure that the information is consistent across the prime and the dual."""
        assert self.prime.k == self.k, "Trivial check k."
        assert self.prime.mesh == self.mesh, "Trivial check mesh."
        assert self.prime.space == self.space, "Trivial check space."
        assert self.prime.orientation == self.orientation, "Trivial check orientation."
        assert self.prime.IS_hybrid is self.IS_hybrid, "prime must be hybrid form, now it is not."
        assert self.prime.IS_volume_form is self.IS_volume_form, "Trivial check IS_volume_form."

    def ___PRIVATE_generate_mass_matrix___(self):
        """For algebra dual forms, this method will only be called once. The result (a mass matrix) will be cached in
        `self._mass_matrix_` because we think it is an essential property for algebra dual forms and it will be used
        for multiple times. Therefore, we cache it somewhere for the dual standard forms.

        :return:  A tuple of two outputs: the mass matrix and the inverse mass matrix.
        """
        MM = self.prime.matrices.mass
        iMM = MM.inv
        return MM, iMM

    @property
    def DO(self):
        """If has too many do methods, we group them in to DO."""
        return self._DO_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def coboundary(self):
        return self._coboundary_







class _3dCSCG_Algebra_DUAL_Standard_Form_DO(FrozenOnly):
    def __init__(self, dsf):
        self._dsf_ = dsf
        self._freeze_self_()



class _3dCSCG_Algebra_DUAL_Standard_Form_Cochain(FrozenOnly):
    """The cochain of algebra dual form is equal to the mass matrix dot the cochain of the prime form.

    dual_cochain = mass_matrix dot prime_cochain.

    """
    def __init__(self, dsf):
        """
        :param dsf: The dual standard form.
        """
        self._dsf_ = dsf
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
            self._local_ = ____3dCSCG_ADSF_Cochain_Local____(self)
            # the local cochain will be renew automatically if the local cochain of the prime form is renewed.
        return self._local_

    @local.setter
    def local(self, local):
        """"""
        assert len(local) == len(self._dsf_.mesh.elements), f"length of local is wrong!"
        prime_local_cochain = dict()
        iM = self._dsf_.inverse_mass_matrix
        for i in self._dsf_.mesh.elements:
            assert i in local, f"mesh element #{i} is missing in local in Core #{rAnk}."
            prime_local_cochain[i] = iM[i] @ local[i]
        self._dsf_.prime.cochain.local = prime_local_cochain

    @property
    def EWC(self):
        """For dual standard forms, this is actually very similar with the local.

        :return:
        """
        return EWC_ColumnVector(self._dsf_.mesh.elements,
                                self.___PRIVATE_local_cochain_link___,
                                cache_key_generator='no_cache')

    def ___PRIVATE_local_cochain_link___(self, i):
        """"""
        MM = self._dsf_.mass_matrix[i]
        LC = self._dsf_.prime.cochain.local[i]
        DLC = csr_matrix(MM @ LC).T
        return DLC


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





class ____3dCSCG_ADSF_Cochain_Local____(FrozenOnly):
    """"""
    def __init__(self, dsf_CO):
        """"""
        self._PC_ = dsf_CO._dsf_.prime.cochain
        self._MM_ = dsf_CO._dsf_.mass_matrix
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






class _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary(FrozenOnly):
    """This is one of the key properties of a standard algebraic dual form. To perform the coboundary, now we need the
    help of another dual trace form.
    """
    def __init__(self, dsf):
        self._dsf_ = dsf
        self._E_ = None
        self._freeze_self_()


    def __call__(self, dual_trace):
        """Compute the dual operator. This must be assisted with a dual trace form which has cochain."""
        assert self._dsf_.k in (1,2,3), f"dual coboundary does not work for 0-dual-standard form."
        assert dual_trace.__class__.__name__ in ('_0_Algebra_DUAL_Trace',
                                                 '_1_Algebra_DUAL_Trace',
                                                 '_2_Algebra_DUAL_Trace'), \
            f"dual trace: {dual_trace} in fact is not a dual trace form."
        assert self._dsf_.mesh == dual_trace.mesh, f"mesh of dual form must match" \
                                                   f"the mesh of the dual trace form."
        assert self._dsf_.space == dual_trace.space, f"space of the dual form does not" \
                                                     f"match that of the dual trace form."

        DFC = self._dsf_.cochain.local # dual standard form cochain
        DTC = dual_trace.cochain.local # dual trace form cochain

        E_T = self.incidence_matrix
        T_T = dual_trace.coboundary.trace_matrix



        if self._dsf_.k == 3: # dual gradient operator: (-E_{div}^T, T_F^T)
            assert dual_trace.k == 2, f"3-dual-form has to match with a 2-dual-trace-form."

            E = - E_T
            T = T_T

            next_prime_form_Path = '_3dCSCG.form.standard._2_form'
            next_prime_form_Name = '_2Form'
            next_dual_form_Path = '_3dCSCG.ADF.standard._2_AD_form'
            next_dual_form_Name = '_2_Algebra_DUAL_Form'

        elif self._dsf_.k == 2: # dual gradient operator: (E_{curl}^T, -T_N^T)
            assert dual_trace.k == 1, f"2-dual-form has to match with a 1-dual-trace-form."

            E = E_T
            T = - T_T

            next_prime_form_Path = '_3dCSCG.form.standard._1_form'
            next_prime_form_Name = '_1Form'
            next_dual_form_Path = '_3dCSCG.ADF.standard._1_AD_form'
            next_dual_form_Name = '_1_Algebra_DUAL_Form'

        elif self._dsf_.k == 1: # dual gradient operator: (-E_{grad}^T, T_N^T)
            assert dual_trace.k == 0, f"1-dual-form has to match with a 0-dual-trace-form."

            E = - E_T
            T = T_T

            next_prime_form_Path = '_3dCSCG.form.standard._0_form'
            next_prime_form_Name = '_0Form'
            next_dual_form_Path = '_3dCSCG.ADF.standard._0_AD_form'
            next_dual_form_Name = '_0_Algebra_DUAL_Form'

        else:
            raise Exception(f"{self._dsf_.k}-dual-form has no coboundary.")



        next_prime_form_class = getattr(import_module(next_prime_form_Path), next_prime_form_Name)
        next_prime_form_Instance = next_prime_form_class(
            self._dsf_.mesh, self._dsf_.space,
            is_hybrid = self._dsf_.prime.IS_hybrid,
            orientation = self._dsf_.orientation,
            numbering_parameters = self._dsf_.prime.numbering._numbering_parameters_,
            name = 'prime-of-dual_operator(' + self._dsf_.standard_properties.name + ')'
        )

        next_dual_form_class = getattr(import_module(next_dual_form_Path), next_dual_form_Name)
        next_dual_form_Instance = next_dual_form_class(
            next_prime_form_Instance,
            next_prime_form_Instance.mesh,
            next_prime_form_Instance.space,
            orientation=next_prime_form_Instance.orientation,
            name='dual_operator(' + self._dsf_.standard_properties.name + ')'
        )

        next_dual_cochain = dict()
        for i in self._dsf_.mesh.elements: # go through all local elements
            next_dual_cochain[i] = E[i] @ DFC[i] + T[i] @ DTC[i]

        next_dual_form_Instance.cochain.local = next_dual_cochain

        return next_dual_form_Instance


    @property
    def incidence_matrix(self):
        """"""
        if self._E_ is None:
            assert self._dsf_.k in (1, 2, 3), \
                f"dual coboundary does not work for 0-dual-standard form."
            if self._dsf_.k == 3:
                E23 = getattr(self._dsf_.space.incidence_matrix, '_2Form').T
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E23, 'constant')
            elif self._dsf_.k == 2:
                E12 = getattr(self._dsf_.space.incidence_matrix, '_1Form').T
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E12, 'constant')
            elif self._dsf_.k == 1:
                E01 = getattr(self._dsf_.space.incidence_matrix, '_0Form').T
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E01, 'constant')
            else:
                raise Exception()

            self._E_ = E

        return self._E_







if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))(
                                        [12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)