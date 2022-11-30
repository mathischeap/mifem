# -*- coding: utf-8 -*-
""""""
import sys
if './' not in sys.path: sys.path.append('./')


from objects.CSCG._3d.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE

from objects.CSCG._3d.ADF.standard.base.do import _3dCSCG_Algebra_DUAL_Standard_Form_DO
from objects.CSCG._3d.ADF.standard.base.cochain.main import _3dCSCG_Algebra_DUAL_Standard_Form_Cochain
from objects.CSCG._3d.ADF.standard.base.coboundary import _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary
from objects.CSCG._3d.ADF.standard.base.whether import _3dCSCG_ADF_SF_Whether
from objects.CSCG._3d.ADF.standard.base.num import _3dCSCG_ADF_SF_NUM
from objects.CSCG._3d.ADF.standard.base.error import _3dCSCG_ADF_SF_Error





class _3dCSCG_Algebra_DUAL_Standard_Form(_3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""

    def __init__(self, ndim, mesh, space, prime, orientation, name):
        """"""
        super().__init__(ndim, mesh, space, prime)
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_form')
        self._DO_ = None
        self._cochain_ = None
        self._coboundary_ = None
        self._whether_ = None
        self._num_ = None
        self._error_ = None


    def ___PreFrozenChecker___(self):
        """This method will be run automatically before we freeze the object. This is very important because
        we have to make sure that the information is consistent across the prime and the dual."""
        assert self.prime.k == self.k, "Trivial check k."
        assert self.prime.mesh == self.mesh, "Trivial check mesh."
        assert self.prime.space == self.space, "Trivial check space."
        assert self.prime.orientation == self.orientation, "Trivial check orientation."
        assert self.prime.whether.hybrid is self.whether.hybrid, "prime must be hybrid form, now it is not."
        assert self.prime.whether.volume_form is self.whether.volume_form, "Trivial check IS_volume_form."

    def ___PRIVATE_generate_mass_matrix___(self):
        """For algebra dual forms, this method will only be called once.

        The result (a mass matrix) will be cached in `self._mass_matrix_` because we think it is an
        essential property for algebra dual forms, and it will be used
        for multiple times. Therefore, we cache it somewhere for the dual standard forms.

        :return:  A tuple of two outputs: the mass matrix and the inverse mass matrix.
        """
        MM = self.prime.matrices.mass
        iMM = MM.inv
        return MM, iMM

    @property
    def orientation(self):
        """An AD standard form can be either inner or outer."""
        return self._orientation_

    @property
    def do(self):
        """If it has too many do methods, we group them in to do."""
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Algebra_DUAL_Standard_Form_DO(self)
        return self._DO_

    @property
    def cochain(self):
        if self._cochain_ is None:
            self._cochain_ = _3dCSCG_Algebra_DUAL_Standard_Form_Cochain(self)
        return self._cochain_

    @property
    def coboundary(self):
        if self._coboundary_ is None:
            self._coboundary_ = _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary(self)
        return self._coboundary_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _3dCSCG_ADF_SF_Whether(self)
        return self._whether_

    @property
    def num(self):
        if self._num_ is None:
            self._num_ = _3dCSCG_ADF_SF_NUM(self)
        return self._num_

    @property
    def error(self):
        if self._error_ is None:
            self._error_ = _3dCSCG_ADF_SF_Error(self)
        return self._error_















if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\base\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))(
                                        [12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)