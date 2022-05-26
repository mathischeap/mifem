# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 1:12 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

# from objects.mpRfT._2d.master import MeshGenerator
from objects.mpRfT._2d.__tests__.Random.mesh import rm
from screws.miscellaneous.mirand import randint
from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mios import remove
from root.read.main import read
from root.save import save


def test_mpRfT2_save_read():
    """"""
    #-------- mesh -----------------------------------------------------------------------
    load = randint(10, 500)
    miprint(f"-sr- [test_mpRfT2_save_read] 2d @ load={load}... ", flush=True)
    mesh2 = rm(load)
    save(mesh2, 'mpRfT2_mesh.mi')
    MESH2 = read('mpRfT2_mesh.mi')
    assert mesh2 == MESH2
    remove('mpRfT2_mesh.mi')

    return 1

if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/__tests__/unittests/save_read.py
    test_mpRfT2_save_read()
