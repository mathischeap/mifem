# -*- coding: utf-8 -*-

from objects.CSCG._2d.master import MeshGenerator as _2dCSCG_MeshGenerator


def ___restore__2dCSCG_Mesh___(parameters, mesh_cache):

    cache_tag = str(parameters)

    if len(mesh_cache) == 0:
        DO_meshing = True
    elif len(mesh_cache) == 2:
        if cache_tag != mesh_cache[0]:
            DO_meshing = True
        else:
            DO_meshing = False
    else:
        raise Exception()

    if DO_meshing:
        assert parameters.pop('type') == '_2dCSCG_Mesh'
        ID = parameters.pop('ID')
        element_layout = parameters.pop('element_layout')
        domain_parameters = parameters.pop('domain_parameters')
        if 'element_distribution_method' in parameters:
            parameters.pop('element_distribution_method')
            EDM = 'debug'  # old distribution methods are all trivial one.
        else:
            EDM = parameters.pop('EDM')
        assert len(parameters) == 0, "make sure all information are used."
        cache_obj = _2dCSCG_MeshGenerator(ID, **domain_parameters)(element_layout, EDM=EDM)
        if len(mesh_cache) == 0:
            mesh_cache.append(cache_tag)
            mesh_cache.append(cache_obj)
        else:
            mesh_cache[0] = cache_tag
            mesh_cache[1] = cache_obj

    else:
        pass

    return mesh_cache[1]
