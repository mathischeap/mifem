# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/2/2022 1:11 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')


from objects.CSCG.tools.gridToVTK import gridToVTK
from objects.CSCG.tools.unstructuredGridToVTK import unstructuredGridToVTK


from objects.CSCG.tools.distribute_local_cochain import distribute_local_cochain




if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/tools/__init__.py
    print(gridToVTK)