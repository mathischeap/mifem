# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 1:49 PM
"""
from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___
from objects.mpRfT._2d.mesh.main import mpRfT2_Mesh


def ___restore__mpRfT2_Mesh___(parameters, mesh_cache):
    """"""
    assert parameters.pop('type') == 'mpRfT2_Mesh'
    cscg_parameters = parameters.pop('cscg')
    cscg = ___restore__2dCSCG_Mesh___(cscg_parameters, mesh_cache)

    dN = parameters.pop('dN')
    refinements = parameters.pop('refinements')

    cscg.___PRIVATE_generate_element_global_numbering___()

    EGN = cscg._element_global_numbering_

    cells_dict = dict()
    for ref in refinements:
        region, indices =  ref.split('|')
        if region in cscg.elements.in_regions:
            ij = indices.split(', ')
            i, j = [int(_) for _ in ij]
            element = EGN[region][i, j]
            if element in cscg.elements:
                for ind in refinements[ref]:
                    if ind == '':
                        rp = str(element)
                    else:
                        rp = str(element) + '-' + ind

                    cells_dict[rp] = refinements[ref][ind]

        else:
            pass

    mesh = mpRfT2_Mesh(cscg, dN, cells_dict)

    return mesh