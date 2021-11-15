
"""Here, we store some most fundamental functions for mifem."""





from SCREWS.miscellaneous import check_filename_mi
from root.config import cHaining, mAster_rank, rAnk, cOmm
import pickle
from _2dCSCG.main import MeshGenerator as _2dCSCG_MeshGenerator
from _2dCSCG.main import SpaceInvoker as _2dCSCG_SpaceInvoker
from _2dCSCG.main import ExactSolutionSelector as _2dCSCG_ExactSolutionSelector
from _2dCSCG.main import FormCaller as _2dCSCG_FormCaller

from _3dCSCG.main import MeshGenerator as _3dCSCG_MeshGenerator
from _3dCSCG.main import SpaceInvoker as _3dCSCG_SpaceInvoker
from _3dCSCG.main import ExactSolutionSelector as _3dCSCG_ExactSolutionSelector
from _3dCSCG.main import FormCaller as _3dCSCG_FormCaller



def save(obj, filename):
    """
    Save a mifem object to file (`.mi` extension).

    :param obj:
    :param filename:
    :return:
    """
    filename = check_filename_mi(filename)
    if isinstance(obj, (list, tuple)):
        _2bs_ = list()
        _sif_ = tuple()
        for obj_i in obj:
            bs = obj_i.___PRIVATE_save___(filename, do_save=False)
            _2bs_.append(bs)

            sif = obj_i.___PRIVATE_saving_info___()
            _sif_ += (sif,)

        _2bs_.append(_sif_)

        cOmm.barrier()
        if rAnk == mAster_rank:
            with open(filename, 'wb') as output:
                pickle.dump(_2bs_, output, pickle.HIGHEST_PROTOCOL)
            output.close()
    else:
        obj.___PRIVATE_save___(filename, do_save=True)






___CACHE_2dCSCG_mesh___ = list() # we only cache the last (ONE) mesh
___CACHE_2dCSCG_space___ = list() # we only cache the last (ONE) space

___CACHE_3dCSCG_mesh___ = list() # we only cache the last (ONE) mesh
___CACHE_3dCSCG_space___ = list() # we only cache the last (ONE) space

def ___PRIVATE_chain_read___(filename):
    if filename[-3:] != '.mi': filename += '.mi'
    with open(filename, 'rb') as inputs:
        obj_dict = pickle.load(inputs)
    inputs.close()
    return obj_dict

def read(filename, read_individuals=None):
    """
    read a mifem object from a file (`.mi` extension).

    :param filename:
    :param read_individuals: For example, if a file contains 5 objects, if we make read_individuals=[1,1,0,1,0], then we
        only read the first, second and fourth objects, for the third and fifth objects, we return None instead.

        If read_individuals is None, we read all objects.
    :return:
    """
    OBJ = cHaining(___PRIVATE_chain_read___, filename)

    if isinstance(OBJ,  list): # multiple saving objects: must be a list.
        if isinstance(OBJ[-1], tuple): # multiple saving objects with saving info: saving info must be in a tuple.
            save_info_tuple = OBJ[-1]
            OBJ = OBJ[:-1]

            for info in save_info_tuple:
                assert info is None or isinstance(info, dict) # saving info for every obj must be None or a dict.

        else:
            pass
    elif isinstance(OBJ, dict): # single saving object: must be a dict.
        OBJ = [OBJ,]
    else:
        raise Exception()

    LEN = len(OBJ)

    if read_individuals is None:
        read_individuals = [1 for _ in range(LEN)]
    else:
        if isinstance(read_individuals, (bool, int, float)):
            read_individuals = [read_individuals,]

        assert isinstance(read_individuals, (list, tuple)), \
            f"read_individuals ({read_individuals.__class__.__name__}) need to be a list or tuple."

        assert len(read_individuals) == LEN, f'file contains {LEN} objects, read_individuals ' \
                                             f'(LEN={len(read_individuals)}) wrong, it must be of len {LEN} '

        for _, ri in enumerate(read_individuals):
            assert isinstance(ri, bool) or ri == 1 or ri == 0, \
                f"read_individuals[{_}]={ri} is wrong, must be True, False, 1 or 0."


    cOmm.barrier()
    objs = tuple()

    for i in range(LEN):
        if read_individuals[i]:
            obj_dict = OBJ[i]
            obj_name = obj_dict.pop('obj')
            obj_para = obj_dict.pop('parameters')
            assert len(obj_dict) == 0, "mi file should store a dict of two keys only."

            if obj_name == '_3dCSCG.mesh.main._3dCSCG_Mesh':
                obj = ___restore__3dCSCG_Mesh___(obj_para, ___CACHE_3dCSCG_mesh___)
            elif obj_name == '_3dCSCG.space.polynomials._3dCSCG_PolynomialSpace':
                obj = ___restore__3dCSCG_Space___(obj_para, ___CACHE_3dCSCG_space___)
            elif obj_name == '_3dCSCG.APP.exact_solutions.main.ExactSolution':
                obj = ___restore__3dCSCG_ExactSolution___(obj_para)
            elif obj_name in ('_3dCSCG.form.standard._0_form._0Form',
                              '_3dCSCG.form.standard._1_form._1Form',
                              '_3dCSCG.form.standard._2_form._2Form',
                              '_3dCSCG.form.standard._3_form._3Form',
                              '_3dCSCG.form.trace._2_trace._2Trace',
                              '_3dCSCG.form.trace._1_trace._1Trace',
                              '_3dCSCG.form.trace._0_trace._0Trace',):
                obj = ___restore__3dCSCG_Form___(obj_para)
            elif obj_name in ('_3dCSCG.ADF.standard._0_AD_form._0_Algebra_DUAL_Form',
                              '_3dCSCG.ADF.standard._1_AD_form._1_Algebra_DUAL_Form',
                              '_3dCSCG.ADF.standard._2_AD_form._2_Algebra_DUAL_Form',
                              '_3dCSCG.ADF.standard._3_AD_form._3_Algebra_DUAL_Form',
                              '_3dCSCG.ADF.trace._0_AD_trace._0_Algebra_DUAL_Trace',
                              '_3dCSCG.ADF.trace._1_AD_trace._1_Algebra_DUAL_Trace',
                              '_3dCSCG.ADF.trace._2_AD_trace._2_Algebra_DUAL_Trace'):
                obj = ___restore__3dCSCG_Algebra_DUAL_Form___(obj_para)
            elif obj_name == '_2dCSCG.mesh.main._2dCSCG_Mesh':
                obj = ___restore__2dCSCG_Mesh___(obj_para, ___CACHE_2dCSCG_mesh___)
            elif obj_name == '_2dCSCG.space.polynomials._2dCSCG_PolynomialSpace':
                obj = ___restore__2dCSCG_Space___(obj_para, ___CACHE_2dCSCG_space___)
            elif obj_name in ('_2dCSCG.form.standard._0_form_inner._0Form_Inner',
                              '_2dCSCG.form.standard._0_form_outer._0Form_Outer',
                              '_2dCSCG.form.standard._1_form_inner._1Form_Inner',
                              '_2dCSCG.form.standard._1_form_outer._1Form_Outer',
                              '_2dCSCG.form.standard._2_form_inner._2Form_Inner',
                              '_2dCSCG.form.standard._2_form_outer._2Form_Outer',
                              '_2dCSCG.form.trace._1_trace_outer._1Trace_Outer',):
                obj = ___restore__2dCSCG_Form___(obj_para)
            else:
                raise Exception(f'Can not restore {obj_name}')

        else:
            obj = None

        objs += (obj,)
        cOmm.barrier()
        OBJ[i] = None # clean memory.

    if len(objs) == 1: objs = objs[0]
    return objs






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
            EDM = 'debug' # old distribution methods are all trivial one.
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

def ___restore__2dCSCG_ExactSolution___(parameters):
    assert parameters.pop('type') == '_2dCSCG_ExactSolution'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    kwargs = parameters.pop('kwargs')
    assert len(parameters) == 0, "make sure all information are used."
    mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, ___CACHE_2dCSCG_mesh___)
    return _2dCSCG_ExactSolutionSelector(mesh)(ID, **kwargs)

def ___restore__2dCSCG_Form___(parameters):
    assert parameters.pop('type') == '_2dCSCG_Form'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    space_parameters = parameters.pop('space_parameters')
    kwargs = parameters.pop('kwargs')

    ___COCHAIN_READ_VERSION___ = -1
    if 'cochain_local' in parameters:
        ___COCHAIN_READ_VERSION___  = 0
        COCHAIN = parameters.pop('cochain_local')
    else:
        ___COCHAIN_READ_VERSION___  = 1
        COCHAIN = parameters.pop('region_wise_cochain_local')

    space = ___restore__2dCSCG_Space___(space_parameters, ___CACHE_2dCSCG_space___)

    TW_func_ES_parameters = parameters.pop('TW_func_ES_parameters')
    TW_current_time = parameters.pop('TW_current_time')

    assert len(parameters) == 0, "make sure all information are used."

    if TW_func_ES_parameters is None:
        # if the func is not from an exact solution, we will not restore the func.
        mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, ___CACHE_2dCSCG_mesh___)
        form = _2dCSCG_FormCaller(mesh, space)(ID, **kwargs)
    else:
        ES_parameters = TW_func_ES_parameters[0]
        ES_variable_name = TW_func_ES_parameters[1]
        ES = ___restore__2dCSCG_ExactSolution___(ES_parameters)


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

    if ___COCHAIN_READ_VERSION___ == 0:
        if COCHAIN != {}:
            cochain_local = dict()
            for i in mesh.elements:
                cochain_local[i] = COCHAIN[i]
            form.cochain.local = cochain_local
    elif ___COCHAIN_READ_VERSION___ == 1:
        form.cochain.___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(COCHAIN)
    else:
        raise Exception()

    return form



def ___restore__3dCSCG_Mesh___(parameters, mesh_cache):
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
        assert parameters.pop('type') == '_3dCSCG_Mesh'
        ID = parameters.pop('ID')
        element_layout = parameters.pop('element_layout')
        domain_parameters = parameters.pop('domain_parameters')
        if 'element_distribution_method' in parameters:
            parameters.pop('element_distribution_method')
            EDM = 'debug' # old distribution methods are all trivial one.
        else:
            EDM = parameters.pop('EDM')
        assert len(parameters) == 0, "make sure all information are used."
        cache_obj = _3dCSCG_MeshGenerator(ID, **domain_parameters)(element_layout, EDM=EDM)
        if len(mesh_cache) == 0:
            mesh_cache.append(cache_tag)
            mesh_cache.append(cache_obj)
        else:
            mesh_cache[0] = cache_tag
            mesh_cache[1] = cache_obj

    else:
        pass

    return mesh_cache[1]

def ___restore__3dCSCG_Space___(parameters, space_cache):
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
        assert parameters.pop('type') == '_3dCSCG_Space'
        ID = parameters.pop('ID')
        inputs = parameters.pop('inputs')
        ndim = parameters.pop('ndim')
        assert len(parameters) == 0, "make sure all information are used."
        cache_obj = _3dCSCG_SpaceInvoker(ID)(inputs, ndim=ndim)
        if len(space_cache) == 0:
            space_cache.append(cache_tag)
            space_cache.append(cache_obj)
        else:
            space_cache[0] = cache_tag
            space_cache[1] = cache_obj

    else:
        pass

    return space_cache[1]

def ___restore__3dCSCG_ExactSolution___(parameters):
    assert parameters.pop('type') == '_3dCSCG_ExactSolution'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    kwargs = parameters.pop('kwargs')
    assert len(parameters) == 0, "make sure all information are used."
    mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, ___CACHE_3dCSCG_mesh___)
    return _3dCSCG_ExactSolutionSelector(mesh)(ID, **kwargs)

def ___restore__3dCSCG_Form___(parameters):
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


    space = ___restore__3dCSCG_Space___(space_parameters, ___CACHE_3dCSCG_space___)

    TW_func_ES_parameters = parameters.pop('TW_func_ES_parameters')
    TW_current_time = parameters.pop('TW_current_time')

    assert len(parameters) == 0, "make sure all information are used."

    if TW_func_ES_parameters is None:

        mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, ___CACHE_3dCSCG_mesh___)
        form = _3dCSCG_FormCaller(mesh, space)(ID, **kwargs)

    else:

        ES_parameters = TW_func_ES_parameters[0]
        ES_variable_name = TW_func_ES_parameters[1]
        ES = ___restore__3dCSCG_ExactSolution___(ES_parameters)




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

def ___restore__3dCSCG_Algebra_DUAL_Form___(parameters):
    """
    For dual forms, we actually directly use the parameters from its prime,
    and the form type is decided by the obj_name.

    """
    assert parameters['type'] == '_3dCSCG_Form'
    ID = parameters['ID']
    prime = ___restore__3dCSCG_Form___(parameters)
    mesh = prime.mesh
    space = prime.space

    if ID == '0-f':
        dID = '0-adf'
    elif ID == '1-f':
        dID = '1-adf'
    elif ID == '2-f':
        dID = '2-adf'
    elif ID == '3-f':
        dID = '3-adf'
    elif ID == '0-t':
        dID = '0-adt'
    elif ID == '1-t':
        dID = '1-adt'
    elif ID == '2-t':
        dID = '2-adt'
    else:
        raise Exception()

    _3FC = _3dCSCG_FormCaller(mesh, space)

    dual = _3FC(dID, prime)

    return dual