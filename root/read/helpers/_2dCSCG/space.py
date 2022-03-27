

from objects.CSCG._2d.master import SpaceInvoker as _2dCSCG_SpaceInvoker

def ___restore__2dCSCG_Space___(parameters, space_cache):
    cache_tag = str(parameters)

    if len(space_cache) == 0:
        DO_meshing = True
    elif len(space_cache) == 2:
        if cache_tag != space_cache[0]:
            DO_meshing = True
        else:
            DO_meshing = False
    else:
        raise Exception()
    if DO_meshing:
        assert parameters.pop('type') == '_2dCSCG_Space'
        ID = parameters.pop('ID')
        inputs = parameters.pop('inputs')
        ndim = parameters.pop('ndim')
        assert len(parameters) == 0, "make sure all information are used."
        cache_obj =  _2dCSCG_SpaceInvoker(ID)(inputs, ndim=ndim)

        if len(space_cache) == 0:
            space_cache.append(cache_tag)
            space_cache.append(cache_obj)
        else:
            space_cache[0] = cache_tag
            space_cache[1] = cache_obj

    else:
        pass

    return space_cache[1]
