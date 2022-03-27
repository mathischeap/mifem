# -*- coding: utf-8 -*-
"""
To test the library, do

mpiexec -n 4 python __tests__/test_all.py

mpiexec -n 4 python __tests__/unittests/main.py

mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/main.py

mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/main.py

"""
import os
absolute_path = os.path.dirname(__file__)
import sys
if absolute_path not in sys.path: sys.path.append(absolute_path)


import objects.CSCG._2d.__init__ as cscg2
import objects.CSCG._3d.__init__ as cscg3


import tools.__init__ as tools
import screws.__init__ as screws


import root.save as save
import root.read.main as read



if __name__ == '__main__':
    print(cscg2)
    print(cscg3)
    print(tools)
    print(screws)
    print(save)
    print(read)

    mesh = cscg2.mesh('rectangle_periodic', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=True)