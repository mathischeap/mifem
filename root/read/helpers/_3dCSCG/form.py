# -*- coding: utf-8 -*-

from objects.CSCG._3d.master import FormCaller as _3dCSCG_FormCaller
from root.read.helpers._3dCSCG.space import ___restore__3dCSCG_Space___
from root.read.helpers._3dCSCG.mesh import ___restore__3dCSCG_Mesh___



def ___restore__3dCSCG_Form___(parameters, mesh_cache, space_cache):
    """"""
    assert parameters.pop('type') == '_3dCSCG_Form'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    space_parameters = parameters.pop('space_parameters')
    kwargs = parameters.pop('kwargs')

    ___COCHAIN_READ_VERSION___ = - 1
    if 'cochain_local' in parameters:
        ___COCHAIN_READ_VERSION___ = 0
        COCHAIN = parameters.pop('cochain_local')
    else:
        ___COCHAIN_READ_VERSION___ = 1
        COCHAIN = parameters.pop('region_wise_cochain_local')


    space = ___restore__3dCSCG_Space___(space_parameters, space_cache)
    mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, mesh_cache)
    form = _3dCSCG_FormCaller(mesh, space)(ID, **kwargs)


    if ___COCHAIN_READ_VERSION___ == 0:
        if COCHAIN != dict():
            cochain_local = dict()
            for i in mesh.elements:
                cochain_local[i] = COCHAIN[i]
            form.cochain.local = cochain_local
    elif ___COCHAIN_READ_VERSION___ == 1:
        form.cochain.___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(COCHAIN)
    else:
        raise Exception()

    return form
