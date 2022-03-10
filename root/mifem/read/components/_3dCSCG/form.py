

from _3dCSCG.main import FormCaller as _3dCSCG_FormCaller
from root.mifem.read.components._3dCSCG.space import ___restore__3dCSCG_Space___
from root.mifem.read.components._3dCSCG.mesh import ___restore__3dCSCG_Mesh___
from root.mifem.read.components._3dCSCG.exact_solution import ___restore__3dCSCG_ExactSolution___



def ___restore__3dCSCG_Form___(parameters, mesh_cache, space_cache):
    """"""
    assert parameters.pop('type') == '_3dCSCG_Form'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    space_parameters = parameters.pop('space_parameters')
    kwargs = parameters.pop('kwargs')

    ___COCHAIN_READ_VERSION___ = - 1
    if 'cochain_local' in parameters:
        ___COCHAIN_READ_VERSION___  = 0
        COCHAIN = parameters.pop('cochain_local')
    else:
        ___COCHAIN_READ_VERSION___  = 1
        COCHAIN = parameters.pop('region_wise_cochain_local')


    space = ___restore__3dCSCG_Space___(space_parameters, space_cache)

    TW_func_ES_parameters = parameters.pop('TW_func_ES_parameters')
    TW_current_time = parameters.pop('TW_current_time')

    assert len(parameters) == 0, "make sure all information are used."

    if TW_func_ES_parameters is None:

        mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, mesh_cache)
        form = _3dCSCG_FormCaller(mesh, space)(ID, **kwargs)

    else:

        ES_parameters = TW_func_ES_parameters[0]
        ES_variable_name = TW_func_ES_parameters[1]
        ES = ___restore__3dCSCG_ExactSolution___(ES_parameters, mesh_cache)




        assert mesh_parameters == ES.mesh.standard_properties.parameters
        mesh = ES.mesh
        # OR use
        # if mesh_parameters == ES.mesh.standard_properties.parameters:
        #     mesh = ES.mesh
        # else:
        #     mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, ___CACHE_3dCSCG_mesh___)





        form = _3dCSCG_FormCaller(mesh, space)(ID, **kwargs)
        form.TW.current_time = TW_current_time
        form.TW.func.___DO_set_func_body_as___(ES, ES_variable_name)
        form.TW.___DO_push_all_to_instant___()

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