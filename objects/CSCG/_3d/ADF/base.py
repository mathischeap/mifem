

""""""

from objects.CSCG.base.ADF.base import CSCG_Algebra_DUAL_FORM_BASE


class _3dCSCG_Algebra_DUAL_FORM_BASE(CSCG_Algebra_DUAL_FORM_BASE):
    """"""

    def __init__(self, ndim, mesh, space):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', "Need a 3dCSCG mesh."
        assert '3dCSCG|structured|space' in space.standard_properties.stamp, "Need a 3dCSCG structured space."
        assert mesh.ndim == space.ndim == 3
        super().__init__(ndim, mesh, space)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_algebra_dual_form')
        assert self.ndim == 3, "CHECK ndim"


    def ___PRIVATE_generate_mass_matrix___(self):
        raise NotImplementedError()