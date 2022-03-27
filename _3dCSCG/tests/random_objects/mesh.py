

from root.config.main import *
import random
from _3dCSCG.master import MeshGenerator



def random_mesh_of_elements_around(elements_num,
                                   exclude_periodic=False,
                                   domain_boundary_distribution_regularities=None,
                                   mesh_boundary_num = None,
                                   mesh_pool=None,
                                   EDM_pool = None,
                                   ):
    """We generate a random mesh of almost ``elements_num`` elements.

    Make sure we edit (1), (2), (3), (4) when new mesh is coded.

    :param elements_num:
    :param exclude_periodic: If we exclude the periodic domains.
    :param mesh_boundary_num: we will select mesh which satisfies this requirement.
    :param domain_boundary_distribution_regularities: We only select from meshes
        that match at least one of `domain_boundary_distribution_regularities`
        (``domain_boundary_distribution_regularity`` \cap mesh.domain.boundaries.distribution_regularities != empty)
    :param mesh_pool: if `mesh_pool` is not None, we will not select mesh from this pool.
    :param EDM_pool: if `EDM_pool` is not None, we will use one of the EDM from this pool. Otherwise, we use
        EDM = None for the mesh generator.
    :return:
    """

    assert isinstance(elements_num, (int, float)) and elements_num >= 1, \
        f"element_num = {elements_num} is wrong, must be a number and >= 1."

    RP = MeshGenerator.___domain_input_random_parameters___()
    Statistic = MeshGenerator.___domain_input_statistic___()


    mesh_name_region_num = dict()
    if exclude_periodic:
        for mesh_id in Statistic:
            if not Statistic[mesh_id]['periodic']:
                mesh_name_region_num[mesh_id] = Statistic[mesh_id]['region num']
    else:
        for mesh_id in Statistic:
            mesh_name_region_num[mesh_id] = Statistic[mesh_id]['region num']



    if mesh_pool is not None:
        if isinstance(mesh_pool, str):
            mesh_pool = [mesh_pool,]
        assert isinstance(mesh_pool, (list, tuple)), f"mesh_pool={mesh_pool} should be a list or a tuple."

        _MESH_POOL_ = dict()

        for mesh in mesh_name_region_num:
            if mesh in mesh_pool:
                _MESH_POOL_[mesh] = mesh_name_region_num[mesh]

        mesh_name_region_num = _MESH_POOL_
    else:
        pass



    mesh_name_boundary_num = dict()
    for mesh_id in Statistic:
        mesh_name_boundary_num[mesh_id] = Statistic[mesh_id]['mesh boundary num']




    if mesh_boundary_num is None:
        pass
    elif len(mesh_boundary_num) >= 3 and mesh_boundary_num[:2] == '>=':
        mesh_boundary_num = int(mesh_boundary_num[2:])
        for mn in mesh_name_region_num:
            if mesh_name_boundary_num[mn] >= mesh_boundary_num:
                pass
            else:
                del mesh_name_region_num[mn]
    else:
        raise NotImplementedError(f"Do not understand mesh_boundary_num={mesh_boundary_num}.")

    mesh_personal_parameters = RP

    assert len(mesh_name_region_num) > 0, f"cannot find a proper mesh."

    #------- check: domain_boundary_distribution_regularities -------------------------------------

    while 1:

        if len(mesh_name_region_num) == 0:
            raise Exception(f"cannot find a mesh with proper domain_boundary_distribution_regularities")

        if rAnk == mAster_rank:
            i = random.randint(0, len(mesh_name_region_num)-1)
            mesh_name = list(mesh_name_region_num.keys())[i]
        else:
            mesh_name = None

        mesh_name = cOmm.bcast(mesh_name, root=mAster_rank)
        if domain_boundary_distribution_regularities is None:
            break

        if isinstance(domain_boundary_distribution_regularities, str):
            domain_boundary_distribution_regularities = \
                [domain_boundary_distribution_regularities,]

        test_mesh = MeshGenerator(mesh_name)([1,1,1])
        dbd_regularities = test_mesh.domain.boundaries.distribution_regularities

        for r in domain_boundary_distribution_regularities:
            if r in dbd_regularities:
                break
        # noinspection PyUnboundLocalVariable
        if r in dbd_regularities: break

        del mesh_name_region_num[mesh_name]



    if mesh_name in mesh_personal_parameters:
        personal_parameters = mesh_personal_parameters[mesh_name]
    else:
        personal_parameters = dict()


    if rAnk == mAster_rank:
        region_num = mesh_name_region_num[mesh_name]
        if elements_num < region_num:
            elements_num = region_num

        elements_num /= region_num
        assert elements_num >= 1

        factor0 = int(elements_num**(1/3))
        if factor0 < 1:
            factor0 = 1

        elements_num /= factor0
        assert elements_num >= 1
        elements_num = np.ceil(elements_num)

        if elements_num == 1:
            factor1 = 1
        else:
            factor1 = random.randint(1,int(elements_num/1.75))

        factor2 = int(elements_num/factor1)
        if factor1 < 1:
            factor1 = 1
        if factor2 < 1:
            factor2 = 1

        #------ (4) EDIT :: special requests for particular meshes ---------------------------------
        if mesh_name == 'crazy_periodic':
            if factor0 == 1: factor0 = 2
            if factor1 == 1: factor1 = 2
            if factor2 == 1: factor2 = 2
        else: # has no special request for the mesh at this moment.
            pass
        #===========================================================================================

        FFF = random.sample((factor0, factor1, factor2), 3)
        element_layout = list()

        a = random.random()

        if a > 0.25: # 75% chance to use non-uniform element_layout
            for f in FFF:
                element_layout.append([random.randint((f+1), 6*(f+1)) for _ in range(f)])
        else: # uniform element_layout
            element_layout = FFF

    else:
        element_layout = None
    element_layout = cOmm.bcast(element_layout, root=mAster_rank)



    if EDM_pool is None:
        EDM = None
    else:
        if rAnk == mAster_rank:

            if isinstance(EDM_pool, str):
                EDM_pool = [EDM_pool,]

            LEN = len(EDM_pool)
            ind = random.randint(0, LEN-1)
            EDM = EDM_pool[ind]
        else:
            EDM = None

        EDM = cOmm.bcast(EDM, root=mAster_rank)


    mesh = MeshGenerator(mesh_name, **personal_parameters)(element_layout, EDM=EDM)
    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        assert mesh.elements.GLOBAL_num == np.prod(FFF) * region_num

    return mesh

