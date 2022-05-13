# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
from objects.nCSCG.rf2._3d.__tests__.Random.mesh import random_mesh_of_elements_around as rm3
from screws.miscellaneous.mirand import randint
from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mios import remove
from root.read.main import read
from root.save import save




def test_nCSCG_RF2_save_read():
    """"""
    load = randint(10, 100)
    miprint(f"-sr- [test_nCSCG_RF2_save_read] 2d @ load={load}... ", flush=True)
    mesh2 = rm2(load)
    save(mesh2, '_2nCSCG_RF2_mesh.mi')
    MESH2 = read('_2nCSCG_RF2_mesh.mi')
    assert mesh2 == MESH2
    remove('_2nCSCG_RF2_mesh.mi')

    load = randint(50, 150)
    miprint(f"     [test_nCSCG_RF2_save_read] 3d @ load={load}... ", flush=True)
    mesh3 = rm3(load)
    save(mesh3, '_3nCSCG_RF2_mesh.mi')
    MESH3 = read('_3nCSCG_RF2_mesh.mi')
    assert mesh3 == MESH3
    remove('_3nCSCG_RF2_mesh.mi')

    return 1




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/base/__tests__/unittests/save_read.py
    test_nCSCG_RF2_save_read()
