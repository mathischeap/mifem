

from root.config.main import cHaining, cOmm

from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___
from root.read.helpers._2dCSCG.space import ___restore__2dCSCG_Space___
from root.read.helpers._2dCSCG.form import ___restore__2dCSCG_Form___

from root.read.helpers._3dCSCG.form import ___restore__3dCSCG_Form___
from root.read.helpers._3dCSCG.mesh import ___restore__3dCSCG_Mesh___
from root.read.helpers._3dCSCG.space import ___restore__3dCSCG_Space___
from root.read.helpers._3dCSCG.ADF import ___restore__3dCSCG_Algebra_DUAL_Form___
from root.read.helpers._3dCSCG.exact_solution import ___restore__3dCSCG_ExactSolution___

from root.read.helpers.chain import chain

___CACHE_2dCSCG_mesh___ = list() # we only cache the last (ONE) mesh
___CACHE_2dCSCG_space___ = list() # we only cache the last (ONE) space

___CACHE_3dCSCG_mesh___ = list() # we only cache the last (ONE) mesh
___CACHE_3dCSCG_space___ = list() # we only cache the last (ONE) space


def read(filename, read_individuals=None):
    """
    read a mifem object from a file (`.mi` extension).

    :param filename:
    :param read_individuals: For example, if a file contains 5 objects, if we make read_individuals=[1,1,0,1,0], then we
        only read the first, second and fourth objects, for the third and fifth objects, we return None instead.

        If read_individuals is None, we read all objects.
    :return:
    """
    OBJ = cHaining(chain, filename)

    if isinstance(OBJ,  list): # multiple saving objects: must be a list.
        if isinstance(OBJ[-1], tuple): # multiple saving objects with saving info: saving info must be in a tuple.
            save_info_tuple = OBJ[-1]
            OBJ = OBJ[:-1]

            for info in save_info_tuple:
                assert info is None or isinstance(info, dict)
                # saving info for every obj must be None or a dict.

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

            if obj_name == '_3dCSCG_Mesh':
                obj = ___restore__3dCSCG_Mesh___(obj_para, ___CACHE_3dCSCG_mesh___)
            elif obj_name == '_3dCSCG_PolynomialSpace':
                obj = ___restore__3dCSCG_Space___(obj_para, ___CACHE_3dCSCG_space___)
            elif obj_name == '_3dCSCG_ExactSolution':
                obj = ___restore__3dCSCG_ExactSolution___(obj_para, ___CACHE_3dCSCG_mesh___)
            elif obj_name in ('_3dCSCG_0Form',
                              '_3dCSCG_1Form',
                              '_3dCSCG_2Form',
                              '_3dCSCG_3Form',
                              '_3dCSCG_2Trace',
                              '_3dCSCG_1Trace',
                              '_3dCSCG_0Trace',
                              '_3dCSCG_0Edge'):
                obj = ___restore__3dCSCG_Form___(obj_para, ___CACHE_3dCSCG_mesh___, ___CACHE_3dCSCG_space___)
            elif obj_name in ('_3dCSCG_S0_ADF',
                              '_3dCSCG_S1_ADF',
                              '_3dCSCG_S2_ADF',
                              '_3dCSCG_S3_ADF',
                              '_3dCSCG_T0_ADF',
                              '_3dCSCG_T1_ADF',
                              '_3dCSCG_T2_ADF'):
                obj = ___restore__3dCSCG_Algebra_DUAL_Form___(obj_para, ___CACHE_3dCSCG_mesh___, ___CACHE_3dCSCG_space___)


            elif obj_name == '_2dCSCG_Mesh':
                obj = ___restore__2dCSCG_Mesh___(obj_para, ___CACHE_2dCSCG_mesh___)
            elif obj_name == '_2dCSCG_PolynomialSpace':
                obj = ___restore__2dCSCG_Space___(obj_para, ___CACHE_2dCSCG_space___)
            elif obj_name in ('_2dCSCG_0Form_Inner',
                              '_2dCSCG_0Form_Outer',
                              '_2dCSCG_1Form_Inner',
                              '_2dCSCG_1Form_Outer',
                              '_2dCSCG_2Form_Inner',
                              '_2dCSCG_2Form_Outer',
                              '_2dCSCG_1Trace_Outer',):
                obj = ___restore__2dCSCG_Form___(obj_para, ___CACHE_2dCSCG_mesh___, ___CACHE_2dCSCG_space___)
            else:
                raise Exception(f'Can not restore {obj_name}')

        else:
            obj = None

        objs += (obj,)
        cOmm.barrier()
        OBJ[i] = None # clean memory.

    if len(objs) == 1: objs = objs[0]
    return objs


