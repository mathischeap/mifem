# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.form.standard._0.base.main import _2nCSCG_RF2_Standard0FormBase


class _2nCSCG_RF2_InnerStandard0Form(_2nCSCG_RF2_Standard0FormBase):
    """"""

    def __init__(self, mesh, hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-0-form'):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
        numbering_parameters
        name
        """
        super(_2nCSCG_RF2_InnerStandard0Form, self).__init__(mesh, hybrid, 'inner', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_standard_inner_0_form')
        self._freeze_self_()





if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/inner/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(50)

    f0 = _2nCSCG_RF2_InnerStandard0Form(mesh)

    from objects.nCSCG.rf2._2d.fields.scalar.main import _2nCSCG_RF2_ScalarField
    import numpy as np
    def p(t, x, y): return np.sin(np.pi * x) * np.cos(np.pi * y) + t
    s = _2nCSCG_RF2_ScalarField(mesh, p)

    f0.TW.func = s
    s.current_time = 0

    f0.discretize()
    f0.visualize()

    print(f0.error.L())

    # mesh.do.update()