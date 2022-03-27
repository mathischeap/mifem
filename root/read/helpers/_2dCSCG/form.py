

from root.read.helpers._2dCSCG.exact_solution import ___restore__2dCSCG_ExactSolution___
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

    TW_func_ES_parameters = parameters.pop('TW_func_ES_parameters')
    TW_current_time = parameters.pop('TW_current_time')

    assert len(parameters) == 0, "make sure all information are used."

    if TW_func_ES_parameters is None:
        # if the func is not from an exact solution, we will not restore the func.
        mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, mesh_cache)
        form = _2dCSCG_FormCaller(mesh, space)(ID, **kwargs)
    else:
        ES_parameters = TW_func_ES_parameters[0]
        ES_variable_name = TW_func_ES_parameters[1]
        ES = ___restore__2dCSCG_ExactSolution___(ES_parameters, mesh_cache)


        assert mesh_parameters == ES.mesh.standard_properties.parameters
        mesh = ES.mesh
        # OR use
        # if mesh_parameters == ES.mesh.standard_properties.parameters:
        #     mesh = ES.mesh
        # else:
        #     mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, ___CACHE_2dCSCG_mesh___)

        form = _2dCSCG_FormCaller(mesh, space)(ID, **kwargs)
        form.TW.current_time = TW_current_time
        form.TW.func.___DO_set_func_body_as___(ES, ES_variable_name)
        form.TW.___DO_push_all_to_instant___()

    form.cochain.___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(COCHAIN)

    return form