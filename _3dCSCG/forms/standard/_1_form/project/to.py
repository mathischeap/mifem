

from screws.frozen import FrozenOnly
from _3dCSCG.forms.standard._0_form.main import _3dCSCG_0Form


class ___3dCSCG_1Form_Project_To___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self,_1sf):
        self._sf_ = _1sf
        self._freeze_self_()

    def standard_2form(self):
        """Project this 1form into a 2form exactly."""
        return self._sf_.special.___PRIVATE_projected_into_2form_exactly___()

    def vector_of_3_standard_0forms(self):
        """project this 1form into a tuple of three 0forms. Each 0form stands for
        a component of the 1form (as a vector).

        Since the 0forms with the same space will be of higher degree, so we do not
        need to make a new space of higher degree. And we of course use the same mesh.
        Thus, both the mesh and space will the same as the those of this 1form."""
        space = self._sf_.space
        mesh = self._sf_.mesh

        f0_x = _3dCSCG_0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
                             orientation=self._sf_.orientation,
                             numbering_parameters=self._sf_.numbering._numbering_parameters_,
                             name='Projected_x_0form_of_'+self._sf_.standard_properties.name)

        f0_y = _3dCSCG_0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
                             orientation=self._sf_.orientation,
                             numbering_parameters=self._sf_.numbering._numbering_parameters_,
                             name='Projected_y_0form_of_'+self._sf_.standard_properties.name)

        f0_z = _3dCSCG_0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
                             orientation=self._sf_.orientation,
                             numbering_parameters=self._sf_.numbering._numbering_parameters_,
                             name='Projected_z_0form_of_'+self._sf_.standard_properties.name)

        _, v = self._sf_.reconstruct(*space.nodes, ravel=True)
        fx_c = dict()
        fy_c = dict()
        fz_c = dict()
        for i in v: # go thorough all local mesh elements
            vx, vy, vz = v[i] # values are actually used as the local cochains of the 0forms

            fx_c[i] = vx
            fy_c[i] = vy
            fz_c[i] = vz

        f0_x.cochain.local = fx_c
        f0_y.cochain.local = fy_c
        f0_z.cochain.local = fz_c

        return f0_x, f0_y, f0_z
