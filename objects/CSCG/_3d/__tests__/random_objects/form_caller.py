"""
Generate random objects for testing purpose.
"""


import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.master import FormCaller
from objects.CSCG._3d.__tests__.random_objects.space import random_space_of_degrees_around
from objects.CSCG._3d.__tests__.random_objects.mesh import random_mesh_of_elements_around
from objects.CSCG._3d.__tests__.random_objects.mesh_and_space import random_mesh_and_space_of_total_load_around




def random_FormCaller_of_total_load_around(*args, **kwargs):
    """A wrapper of `random_mesh_and_space_of_total_load_around` and we use the outputs to make a
    3D FormCaller instance."""
    mesh, space = random_mesh_and_space_of_total_load_around(*args, **kwargs)
    return FormCaller(mesh, space)





if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\TESTS\random_objects.py
    # random_mesh_of_elements_around(1)
    M = random_mesh_of_elements_around
    S = random_space_of_degrees_around

    random_mesh_and_space_of_total_load_around(100)
