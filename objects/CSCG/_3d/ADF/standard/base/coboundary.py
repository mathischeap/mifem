# -*- coding: utf-8 -*-
from importlib import import_module
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix


class _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary(FrozenOnly):
    """This is one of the key properties of a algebraic dual form.

    To perform the coboundary, now we need the help of another dual trace form.
    """
    def __init__(self, dsf):
        self._dsf_ = dsf
        self._E_ = None
        self._freeze_self_()

    def __call__(self, dual_trace):
        """Compute the dual operator. This must be assisted with a dual trace form which has cochain."""
        assert self._dsf_.k in (1, 2, 3), f"dual coboundary does not work for 0-dual-standard form."
        assert dual_trace.__class__.__name__ in ('_3dCSCG_T0_ADF',
                                                 '_3dCSCG_T1_ADF',
                                                 '_3dCSCG_T2_ADF'), \
            f"dual trace: {dual_trace} in fact is not a dual trace form."
        assert self._dsf_.mesh == dual_trace.mesh, f"mesh of dual form must match" \
                                                   f"the mesh of the dual trace form."
        assert self._dsf_.space == dual_trace.space, f"space of the dual form does not" \
                                                     f"match that of the dual trace form."

        DFC = self._dsf_.cochain.local  # dual standard form cochain
        DTC = dual_trace.cochain.local  # dual trace form cochain

        E_T = self.incidence_matrix.T
        T_T = dual_trace.coboundary.trace_matrix.T

        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

        if self._dsf_.k == 3:  # dual gradient operator: (-E_{div}^T, T_F^T)
            assert dual_trace.k == 2, f"3-dual-form has to match with a 2-dual-trace-form."

            E = - E_T
            T = T_T

            next_prime_form_Path = base_path + 'forms.standard._2s.main'
            next_prime_form_Name = '_3dCSCG_2Form'
            next_dual_form_Path = base_path + 'ADF.standard._2s.main'
            next_dual_form_Name = '_3dCSCG_S2_ADF'

        elif self._dsf_.k == 2:  # dual gradient operator: (E_{curl}^T, -T_N^T)
            assert dual_trace.k == 1, f"2-dual-form has to match with a 1-dual-trace-form."

            E = E_T
            T = - T_T

            next_prime_form_Path = base_path + 'forms.standard._1s.main'
            next_prime_form_Name = '_3dCSCG_1Form'
            next_dual_form_Path = base_path + 'ADF.standard._1s.main'
            next_dual_form_Name = '_3dCSCG_S1_ADF'

        elif self._dsf_.k == 1:  # dual gradient operator: (-E_{grad}^T, T_N^T)
            assert dual_trace.k == 0, f"1-dual-form has to match with a 0-dual-trace-form."

            E = - E_T
            T = T_T

            next_prime_form_Path = base_path + 'forms.standard._0s.main'
            next_prime_form_Name = '_3dCSCG_0Form'
            next_dual_form_Path = base_path + 'ADF.standard._0s.main'
            next_dual_form_Name = '_3dCSCG_S0_ADF'

        else:
            raise Exception(f"{self._dsf_.k}-dual-form has no coboundary.")

        next_prime_form_class = getattr(import_module(next_prime_form_Path), next_prime_form_Name)
        next_prime_form_Instance = next_prime_form_class(
            self._dsf_.mesh, self._dsf_.space,
            hybrid=self._dsf_.prime.whether.hybrid,
            orientation=self._dsf_.orientation,
            numbering_parameters=self._dsf_.prime.numbering._numbering_parameters_,
            name='prime-of-dual_operator(' + self._dsf_.standard_properties.name + ')'
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
        for i in self._dsf_.mesh.elements:  # go through all local elements
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
                E32 = getattr(self._dsf_.space.incidence_matrix, '_3dCSCG_2Form')
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E32, 'constant')
            elif self._dsf_.k == 2:
                E21 = getattr(self._dsf_.space.incidence_matrix, '_3dCSCG_1Form')
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E21, 'constant')
            elif self._dsf_.k == 1:
                E10 = getattr(self._dsf_.space.incidence_matrix, '_3dCSCG_0Form')
                E = EWC_SparseMatrix(self._dsf_.mesh.elements, E10, 'constant')
            else:
                raise Exception()

            self._E_ = E

        return self._E_
