# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 6:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')


from root.read.helpers._3dCSCG.mesh import ___restore__3dCSCG_Mesh___
from objects.nCSCG.rf2._3d.mesh.main import _3nCSCG_RF2_Mesh


def ___restore__3nCSCG_RF2_Mesh___(parameters, mesh_cache):
    """"""
    assert parameters.pop('type') == '_3nCSCG_RF2_Mesh'
    cscg_parameters = parameters.pop('cscg')
    cscg = ___restore__3dCSCG_Mesh___(cscg_parameters, mesh_cache)

    mesh = _3nCSCG_RF2_Mesh(cscg)

    refinements = parameters['refinements']

    for ref in refinements:
        lv0 = ref[0]
        if lv0 in cscg.elements:
            _ = mesh(ref, dynamic=True)

    return mesh


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
