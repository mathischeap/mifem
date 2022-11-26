# -*- coding: utf-8 -*-

from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___
from root.read.helpers._2dCSCG.space import ___restore__2dCSCG_Space___
from objects.CSCG._2d.master import FormCaller as _2dCSCG_FormCaller

def ___restore__2dCSCG_Form___(parameters, mesh_cache, space_cache):
    assert parameters.pop('type') == '_2dCSCG_Form'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    space_parameters = parameters.pop('space_parameters')
    kwargs = parameters.pop('kwargs')

    COCHAIN = parameters.pop('region_wise_cochain_local')

    space = ___restore__2dCSCG_Space___(space_parameters, space_cache)
    mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, mesh_cache)
    form = _2dCSCG_FormCaller(mesh, space)(ID, **kwargs)


    form.cochain.___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(COCHAIN)

    return form