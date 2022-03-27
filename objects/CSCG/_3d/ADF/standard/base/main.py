
"""

"""
import sys
if './' not in sys.path: sys.path.append('/')


from objects.CSCG._3d.ADF.base import _3dCSCG_Algebra_DUAL_FORM_BASE

from objects.CSCG.base.ADF.standard.main import CSCG_Algebra_DUAL_Standard_Form

from objects.CSCG._3d.ADF.standard.base.do import _3dCSCG_Algebra_DUAL_Standard_Form_DO
from objects.CSCG._3d.ADF.standard.base.cochain.main import _3dCSCG_Algebra_DUAL_Standard_Form_Cochain
from objects.CSCG._3d.ADF.standard.base.coboundary import _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary



class _3dCSCG_Algebra_DUAL_Standard_Form(CSCG_Algebra_DUAL_Standard_Form, _3dCSCG_Algebra_DUAL_FORM_BASE):
    """"""

    def __init__(self, ndim, mesh, space, orientation, name):
        """"""
        super().__init__(ndim, mesh, space)
        super().___init___()
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_form')
        self._DO_ = _3dCSCG_Algebra_DUAL_Standard_Form_DO(self)
        self._cochain_ = _3dCSCG_Algebra_DUAL_Standard_Form_Cochain(self)
        self._coboundary_ = _3dCSCG_Algebra_DUAL_Standard_Form_Coboundary(self)

    def ___PRIVATE_reset_cache___(self):
        """"""

    def ___PreFrozenChecker___(self):
        """This method will be run automatically before we freeze the object. This is very important because
        we have to make sure that the information is consistent across the prime and the dual."""
        assert self.prime.k == self.k, "Trivial check k."
        assert self.prime.mesh == self.mesh, "Trivial check mesh."
        assert self.prime.space == self.space, "Trivial check space."
        assert self.prime.orientation == self.orientation, "Trivial check orientation."
        assert self.prime.IS.hybrid is self.IS.hybrid, "prime must be hybrid form, now it is not."
        assert self.prime.IS.volume_form is self.IS.volume_form, "Trivial check IS_volume_form."

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
        """If has too many do methods, we group them in to do."""
        return self._DO_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def coboundary(self):
        return self._coboundary_
















if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\base\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))(
                                        [12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)