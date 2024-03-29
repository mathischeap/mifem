# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/02 11:00 PM
"""
import numpy as np
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.miscellaneous.miprint import miprint
from components.miscellaneous.mirand import randint
from components.miscellaneous.mios import remove
from objects.CSCG.tools.__init__ import unstructuredGridToVTK
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around as rf
from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around as rf2
from tests.objects.CSCG._3d.randObj.field import random_scalar, random_vector
from tests.objects.CSCG._2d.randObj.field import random_scalar as rs2
from tests.objects.CSCG._2d.randObj.field import random_vector as rv2


def TEST_save_CSCG_objects_to_unstructured_VTK_file():
    miprint("VTK [TEST_save_CSCG_objects_to_unstructured_VTK_file] ....", flush=True)

    FC = rf(randint(99, 199))

    scalar = random_scalar(FC.mesh)
    velocity = random_vector(FC.mesh)

    f0 = FC('0-f', hybrid=False, name='pressure')
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    f1 = FC('1-f', hybrid=False, name='vorticity')
    f1.CF = velocity
    f1.CF.current_time = 0
    f1.discretize()
    f2 = FC('2-f', hybrid=False, name='velocity')
    f2.CF = velocity
    f2.CF.current_time = 0
    f2.discretize()
    f3 = FC('3-f', hybrid=False, name='total pressure')
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()

    grid = [np.linspace(-1, 1, 10), np.linspace(-1, 1, 11), np.linspace(-1, 1, 9)]
    filename = '__unstructuredGridToVTK_test__'
    unstructuredGridToVTK(grid, [f0, f1, f2, f3], filename)
    remove(filename + '.vtu')

    FC = rf2(randint(99, 199))
    scalar = rs2(FC.mesh)
    velocity = rv2(FC.mesh)
    f0 = FC('0-f-o', hybrid=False, name='pressure')
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    f1i = FC('1-f-i', hybrid=False, name='Velocity')
    f1i.CF = velocity
    f1i.CF.current_time = 0
    f1i.discretize()
    f1o = FC('1-f-o', hybrid=False, name='vorticity')
    f1o.CF = velocity
    f1o.CF.current_time = 0
    f1o.discretize()
    f2 = FC('2-f-o', hybrid=False, name='total pressure')
    f2.CF = scalar
    f2.CF.current_time = 0
    f2.discretize()

    grid = [np.linspace(-1, 1, 20), np.linspace(-1, 1, 19)]
    filename = '__unstructuredGridToVTK_2dtest__'
    unstructuredGridToVTK(grid, [f0, f1o, f1i, f2], filename)
    remove(filename + '.vtu')

    return 1


if __name__ == "__main__":
    # mpiexec -n 4 python tests/tools/VTK/unstructuredGridToVTK.py
    TEST_save_CSCG_objects_to_unstructured_VTK_file()
