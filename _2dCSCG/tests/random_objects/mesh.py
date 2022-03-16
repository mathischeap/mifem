
from root.config.main import *
import random

from _2dCSCG.main import MeshGenerator


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


    #------(1) EDIT the dict -> mesh_name_region_num :: we will select from this dict ----------------
    if exclude_periodic:
        mesh_name_region_num = {  # add new meshes to this pool.
            # mesh name: how many regions
            'chp1': 4,
            'chp2': 8,
            'crazy': 1,
            'quadrangle': 1,
            'bcr': 1,
            'cic': 6,
            'rectangle': 'unknown'
        }
    else:
        mesh_name_region_num = { # add new meshes to this pool.
            # mesh name: how many regions
            'chp1': 4,
            'chp2': 8,
            'crazy': 1,
            'quadrangle': 1,
            'crazy_periodic': 1,
            'bcr': 1,
            'cic': 6,
            'rectangle': 'unknown'
        }

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


    #---------(2) EDIT the dict -> mesh_name_boundary_num :: parse mesh_boundary_num -----------------
    mesh_name_boundary_num = {
            'chp1': 5,
            'chp2': 5,
            'crazy': 4,
            'quadrangle': 4,
            'crazy_periodic': 0,
            'bcr': 4,
            'cic': 5,
            'rectangle': 4
    }

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

    # -------(3) EDIT the dict -> mesh_personal_parameters :: random personal parameters ------------
    if rAnk == mAster_rank:
        mesh_personal_parameters = { # edit below when new mesh is added to above pool.
            'crazy':{'c': random.randint(0,3)*random.random()/10,
                     'bounds': [(-random.random(), random.random()+0.5),
                                (-random.random(), random.random()+0.5)]
                      },
            'crazy_periodic':{'c': random.randint(0,3)*random.random()/10,
                              'bounds': [(-random.random(), random.random()+0.5),
                                         (-random.random(), random.random()+0.5)]
                              },
            'rectangle': {'p_UL': (random.uniform(-1,1), random.uniform(-1,1)),
                          'width':random.uniform(1,3),
                          'length':random.uniform(2,4),
                          "region_layout": (random.randint(1,3), random.randint(1,3))
                          }
        }
    else:
        mesh_personal_parameters = None

    mesh_personal_parameters = cOmm.bcast(mesh_personal_parameters, root=mAster_rank)

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

        test_mesh = MeshGenerator(mesh_name)([1,1])
        dbd_regularities = test_mesh.domain.boundaries.distribution_regularities

        for r in domain_boundary_distribution_regularities:
            if r in dbd_regularities:
                break
        # noinspection PyUnboundLocalVariable
        if r in dbd_regularities: break

        del mesh_name_region_num[mesh_name]

    #===========================================================

    if mesh_name in mesh_personal_parameters:
        personal_parameters = mesh_personal_parameters[mesh_name]
    else:
        personal_parameters = dict()

    test_mesh = MeshGenerator(mesh_name, **personal_parameters)([1,1])

    #===========================================================

    if rAnk == mAster_rank:
        region_num = mesh_name_region_num[mesh_name]

        if region_num == 'unknown':
            region_num = test_mesh.domain.regions.num
        elif isinstance(region_num, int):
            pass
        else:
            raise Exception(f"region_num={region_num} wrong, can only be positive integer or 'unknown")

        if elements_num < region_num:
            elements_num = region_num

        elements_num /= region_num
        assert elements_num >= 1

        factor0 = int(elements_num**(1/2))
        if factor0 < 1:
            factor0 = 1

        elements_num /= factor0
        assert elements_num >= 1
        elements_num = np.ceil(elements_num)

        if elements_num == 1:
            factor1 = 1
        else:
            factor1 = random.randint(1,elements_num)


        #------ (4) EDIT :: special requests for particular meshes ---------------------------------
        if mesh_name == 'crazy_periodic':
            # factor0 and factor1 are elements along each direction.
            if factor0 == 1: factor0 = 2
            if factor1 == 1: factor1 = 2
        else: # has no special request for the mesh at this moment.
            pass
        #===========================================================================================

        FFF = random.sample((factor0, factor1), 2)

        element_layout = list()
        for f in FFF:
            element_layout.append([random.randint((f+1), 6*(f+1)) for _ in range(f)])

    else:
        element_layout = None
    element_layout = cOmm.bcast(element_layout, root=mAster_rank)

    #-----------------------------------------------------------------------------------------------

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

    # now make the mesh ----------------------------------------------------------------------------
    mesh = MeshGenerator(mesh_name, **personal_parameters)(element_layout, EDM=EDM)

    if rAnk == mAster_rank:
        # noinspection PyUnboundLocalVariable
        assert mesh.elements.GLOBAL_num == np.prod(FFF) * region_num

    return mesh
