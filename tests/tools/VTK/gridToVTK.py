# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/02 10:28 PM
"""
import numpy as np
import sys

if './' not in sys.path: sys.path.append('./')
from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mirand import randint
from screws.miscellaneous.mios import cleandir, rmdir
from objects.CSCG.tools.__init__ import gridToVTK
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around as rf3
from tests.objects.CSCG._3d.randObj.field import random_scalar, random_vector


from tests.objects.CSCG._2d.randObj.form_caller import random_FormCaller_of_total_load_around as rf2
from tests.objects.CSCG._2d.randObj.field import random_scalar as rs2
from tests.objects.CSCG._2d.randObj.field import random_vector as rv2


def TEST_save_CSCG_objects_to_structured_VTK_file():
    miprint("VTK [TEST_save_CSCG_objects_to_structured_VTK_file] ....", flush=True)

    FC = rf3(randint(99, 199))

    scalar = random_scalar(FC.mesh)
    velocity = random_vector(FC.mesh)

    f0 = FC('0-f', is_hybrid=False, name='pressure')
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    f1 = FC('1-f', is_hybrid=False, name='vorticity')
    f1.CF = velocity
    f1.CF.current_time = 0
    f1.discretize()
    f2 = FC('2-f', is_hybrid=False, name='velocity')
    f2.CF = velocity
    f2.CF.current_time = 0
    f2.discretize()
    f3 = FC('3-f', is_hybrid=False, name='total pressure')
    f3.CF = scalar
    f3.CF.current_time = 0
    f3.discretize()

    grid = [np.linspace(-1,1,11), np.linspace(-1,1,8), np.linspace(-1,1,9)]

    path = '__gridToVTK_test__'
    gridToVTK(grid, [f0, f1, f2, f3], path)
    cleandir(path)
    rmdir(path)

    # ----------- 2d test below ---------------------------------------------------------
    FC2 = rf2(randint(99, 199))
    scalar = rs2(FC2.mesh)
    velocity = rv2(FC2.mesh)

    f0 = FC2('0-f-o', is_hybrid=False, name='pressure')
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()
    f1i = FC2('1-f-i', is_hybrid=False, name='Velocity')
    f1i.CF = velocity
    f1i.CF.current_time = 0
    f1i.discretize()
    f1o = FC2('1-f-o', is_hybrid=False, name='vorticity')
    f1o.CF = velocity
    f1o.CF.current_time = 0
    f1o.discretize()
    f2 = FC2('2-f-o', is_hybrid=False, name='total pressure')
    f2.CF = scalar
    f2.CF.current_time = 0
    f2.discretize()

    grid = [np.linspace(-1,1,23), np.linspace(-1,1,17)]
    path = '__gridToVTK_2dtest__'
    gridToVTK(grid, [f0, f1i, f1o, f2], path)
    cleandir(path)
    rmdir(path)

    return 1

if __name__ == "__main__":
    # mpiexec -n 4 python __tests__/unittests/VTK/gridToVTK.py
    TEST_save_CSCG_objects_to_structured_VTK_file()
