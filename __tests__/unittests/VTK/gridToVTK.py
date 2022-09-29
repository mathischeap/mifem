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
from objects.CSCG._3d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around as rf3
from objects.CSCG._3d.__tests__.Random.field import random_scalar, random_vector


from objects.CSCG._2d.__tests__.Random.form_caller import random_FormCaller_of_total_load_around as rf2
from objects.CSCG._2d.__tests__.Random.field import random_scalar as rs2
from objects.CSCG._2d.__tests__.Random.field import random_vector as rv2


def TEST_save_CSCG_objects_to_structured_VTK_file():
    miprint("VTK [TEST_save_CSCG_objects_to_structured_VTK_file] ....", flush=True)

    FC = rf3(randint(99, 199))

    scalar = random_scalar(FC.mesh)
    velocity = random_vector(FC.mesh)

    f0 = FC('0-f', is_hybrid=False, name='pressure')
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.___DO_push_all_to_instant___()
    f0.discretize()
    f1 = FC('1-f', is_hybrid=False, name='vorticity')
    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()
    f2 = FC('2-f', is_hybrid=False, name='velocity')
    f2.TW.func.do.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()
    f3 = FC('3-f', is_hybrid=False, name='total pressure')
    f3.TW.func.do.set_func_body_as(scalar)
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
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
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.___DO_push_all_to_instant___()
    f0.discretize()
    f1i = FC2('1-f-i', is_hybrid=False, name='Velocity')
    f1i.TW.func.do.set_func_body_as(velocity)
    f1i.TW.current_time = 0
    f1i.TW.___DO_push_all_to_instant___()
    f1i.discretize()
    f1o = FC2('1-f-o', is_hybrid=False, name='vorticity')
    f1o.TW.func.do.set_func_body_as(velocity)
    f1o.TW.current_time = 0
    f1o.TW.___DO_push_all_to_instant___()
    f1o.discretize()
    f2 = FC2('2-f-o', is_hybrid=False, name='total pressure')
    f2.TW.func.do.set_func_body_as(scalar)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
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
