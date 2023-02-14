# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard._0s.main import _3dCSCG_0Form
from objects.CSCG._3d.forms.standard._2s.main import _3dCSCG_2Form


class ___3dCSCG_1Form_Project_To___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._freeze_self_()

    def standard_2form(self):
        """Project this 1form into a 2form exactly."""
        space = self._sf_.space
        mesh = self._sf_.mesh

        sp = space.p
        op = list()
        for i in sp:
            op.append(i+1)

        SPACE = space.__class__(op, None)

        f2 = _3dCSCG_2Form(
            mesh, SPACE,
            hybrid=self._sf_.whether.hybrid,
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

    def vector_of_3_standard_0forms(self):
        """project this 1form into a tuple of three 0forms. Each 0form stands for
        a component of the 1form (as a vector).

        Since the 0forms with the same space will be of higher degree, so we do not
        need to make a new space of higher degree. And we of course use the same mesh.
        Thus, both the mesh and space will the same as the those of this 1form."""
        space = self._sf_.space
        mesh = self._sf_.mesh

        f0_x = _3dCSCG_0Form(
            mesh, space, hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_x_0form_of_'+self._sf_.standard_properties.name
        )

        f0_y = _3dCSCG_0Form(
            mesh, space, hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_y_0form_of_'+self._sf_.standard_properties.name
        )

        f0_z = _3dCSCG_0Form(
            mesh, space, hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_z_0form_of_'+self._sf_.standard_properties.name
        )

        _, v = self._sf_.reconstruct(*space.nodes, ravel=True)
        fx_c = dict()
        fy_c = dict()
        fz_c = dict()
        for i in v:  # go thorough all local mesh elements
            vx, vy, vz = v[i]  # values are actually used as the local cochains of the 0forms

            fx_c[i] = vx
            fy_c[i] = vy
            fz_c[i] = vz

        f0_x.cochain.local = fx_c
        f0_y.cochain.local = fy_c
        f0_z.cochain.local = fz_c

        return f0_x, f0_y, f0_z
