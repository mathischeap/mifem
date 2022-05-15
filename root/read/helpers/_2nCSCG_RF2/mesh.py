

from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___
from objects.nCSCG.rf2._2d.mesh.main import _2nCSCG_RF2_Mesh


def ___restore__2nCSCG_RF2_Mesh___(parameters, mesh_cache):
    """"""
    assert parameters.pop('type') == '_2nCSCG_RF2_Mesh'
    cscg_parameters = parameters.pop('cscg')
    cscg = ___restore__2dCSCG_Mesh___(cscg_parameters, mesh_cache)

    space_parameters = parameters.pop('space')

    dN, space_type, space_kwargs = space_parameters
    mesh = _2nCSCG_RF2_Mesh(cscg, dN, space_type=space_type, **space_kwargs)

    refinements = parameters.pop('refinements')

    cscg.___PRIVATE_generate_element_global_numbering___()

    EGN = cscg._element_global_numbering_

    mesh.do.unlock()
    for ref in refinements:
        region, indices =  ref.split('|')
        if region in cscg.elements.in_regions:
            ij = indices.split(', ')
            i, j = [int(_) for _ in ij]
            element = EGN[region][i, j]
            if element in cscg.elements:
                for ind in refinements[ref]:
                    IND = (element,) + ind
                    mesh(IND, dynamic=True)
        else:
            pass

    mesh.do.update() # update the mesh.
    cscg.___element_global_numbering___ = None

    return mesh