# -*- coding: utf-8 -*-

import sys
if './' not in sys.path:
    sys.path.append('/')

from objects.CSCG._2d.master import FormCaller

from tests.objects.CSCG._2d.randObj.mesh_and_space import random_mesh_and_space_of_total_load_around, \
    random_mesh_of_elements_around, random_space_of_degrees_around


def random_FormCaller_of_total_load_around(*args, **kwargs):
    """A wrapper of `random_mesh_and_space_of_total_load_around` and we use the outputs to make a
    3D FormCaller instance."""
    mesh, space = random_mesh_and_space_of_total_load_around(*args, **kwargs)
    return FormCaller(mesh, space)


if __name__ == '__main__':
    # mpiexec -n 4 python _2dCSCG\tests\random_objects.py
    # random_mesh_of_elements_around(1)
    S = random_space_of_degrees_around
    M = random_mesh_of_elements_around

    FC = random_FormCaller_of_total_load_around(100)
    mesh = FC._mesh_
    space = FC._space_
    print(mesh, mesh.elements.layout, space.p)
