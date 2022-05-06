

from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___
from objects.nCSCG.rf2._2d.mesh.main import _2nCSCG_RF2_Mesh


def ___restore__2nCSCG_RF2_Mesh___(parameters, mesh_cache):
    """"""
    assert parameters.pop('type') == '_2nCSCG_RF2_Mesh'
    cscg_parameters = parameters.pop('cscg')
    cscg = ___restore__2dCSCG_Mesh___(cscg_parameters, mesh_cache)

    mesh = _2nCSCG_RF2_Mesh(cscg)

    refinements = parameters['refinements']

    for ref in refinements:
        lv0 = ref[0]
        if lv0 in cscg.elements:
            _ = mesh(ref, dynamic=True)

    return mesh