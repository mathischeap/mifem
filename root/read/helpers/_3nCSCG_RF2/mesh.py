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

    refinements = parameters.pop('refinements')

    cscg.___PRIVATE_generate_element_global_numbering___()

    EGN = cscg._element_global_numbering_

    mesh.do.unlock()
    for ref in refinements:
        region, indices =  ref.split('|')
        if region in cscg.elements.in_regions:
            ijk = indices.split(', ')
            i, j, k = [int(_) for _ in ijk]
            element = EGN[region][i, j, k]
            if element in cscg.elements:
                for ind in refinements[ref]:
                    IND = (element,) + ind
                    mesh(IND, dynamic=True)
        else:
            pass

    mesh.do.update()
    cscg.___element_global_numbering___ = None

    return mesh


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
