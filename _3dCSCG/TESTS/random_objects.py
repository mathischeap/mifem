"""

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
import random
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector






def random_3D_mesh_of_elements_around(elements_num,
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
            'crazy': 1,
            'bridge_arch_cracked': 4,
        }
    else:
        mesh_name_region_num = { # add new meshes to this pool.
            # mesh name: how many regions
            'crazy': 1,
            'crazy_periodic': 1,
            'bridge_arch_cracked': 4,
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
            'crazy': 6,
            'crazy_periodic': 0, # crazy_periodic mesh has no mesh boundary.
            'bridge_arch_cracked': 8,
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
                                (-random.random(), random.random()+0.5),
                                (-random.random(), random.random()+0.5)]
                      },

            'crazy_periodic': {'c': random.randint(0,3)*random.random()/10,
                               'bounds': [(-random.random(), random.random()+0.5),
                                          (-random.random(), random.random()+0.5),
                                          (-random.random(), random.random()+0.5)]
                               },
        }
    else:
        mesh_personal_parameters = None

    mesh_personal_parameters = cOmm.bcast(mesh_personal_parameters, root=mAster_rank)


    # ==============================================================================================

    while 1:

        if len(mesh_name_region_num) == 0:
            raise Exception(f"cannot find a proper mesh.")

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
        for f in FFF:
            element_layout.append([random.randint((f+1), 6*(f+1)) for _ in range(f)])

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



def random_3D_space_of_overall_degrees_around(overall_degrees):
    """

    :param overall_degrees: p0 * p1 * p2
    :return:
    """
    assert isinstance(overall_degrees, (int, float)) and overall_degrees >= 1, \
        f"overall_degrees = {overall_degrees} is wrong, must be a number and >= 1."
    if rAnk == mAster_rank:

        space_pool = ('polynomials',)
        space_name = random.sample(space_pool, 1)[0]

        overall_degrees = int(overall_degrees)

        if overall_degrees == 1:
            p0 = p1 = p2 = 1
        elif overall_degrees == 2:
            p0 = p1 = 1
            p2 = 2
        elif overall_degrees == 3:
            p0 = p1 = 1
            p2 = 3
        elif overall_degrees == 4:
            p0 = p1 = 2
            p2 = 1
        elif overall_degrees == 5:
            p0 = p1 = 2
            p2 = 1
        elif overall_degrees == 6:
            p0 = 3
            p1 = 2
            p2 = 1
        else:
            p0 = random.randint(1, int(overall_degrees/3))
            p1p2 = np.ceil(overall_degrees / p0)
            p1 = random.randint(1, int(p1p2/2))
            p2 = int(p1p2/p1)

        if p0 < 1: p0 = 1
        if p1 < 1: p1 = 1
        if p2 < 1: p2 = 1

        PPP = random.sample((p0, p1, p2), 3)

        DIS = list()
        for p in PPP:

            dis = random.sample(range(p+1, 5*(p+1)), p+1)
            dis.sort()
            dis = np.array(dis)
            dis = dis - np.min(dis)
            dis = dis / np.max(dis)
            dis = dis * 2 - 1
            dis[0] = -1
            dis[-1] = 1
            DIS.append([dis,])
    else:
        space_name = None
        DIS = None
    space_name, DIS = cOmm.bcast([space_name, DIS], root=mAster_rank)

    space = SpaceInvoker(space_name)(DIS)

    return space



def random_3D_mesh_and_space_of_total_load_around(total_load, **kwargs):
    """

    :param total_load: K * P, we will try to make sure K > P.
        K: elements_num
        P: overall_degrees
    :param kwargs:
    :return:
    """
    assert isinstance(total_load, (int, float)) and total_load >= 1, \
        f"total_load = {total_load} is wrong, must be a number and >= 1."

    if rAnk == mAster_rank:
        if total_load == 1:
            K, P = 1, 1
        elif total_load == 2:
            K, P = 2, 1
        elif total_load == 3:
            K, P = 3, 1
        elif total_load == 4:
            K, P = [[2,2],[4,1]][random.randint(0,1)]
        elif total_load == 5:
            K, P = [[3,2], [5,1]][random.randint(0,1)]
        elif total_load == 6:
            K, P = [[3,2], [6,1], [4, 2]][random.randint(0,2)]
        elif total_load == 7:
            K, P = [[3,2], [7,1], [4, 2]][random.randint(0,2)]
        elif total_load == 8:
            K, P = [[3,2], [7,1], [4, 2], [3,3]][random.randint(0,3)]
        else:
            f0 = random.randint(int(total_load**0.5), int(3*total_load**0.5))
            f1 = int(total_load/f0)
            f0, f1 = random.sample((f0, f1), 2)

            if f0 > f1:
                K = f0
                P = f1
            else:
                P = f0
                K = f1

            if K < 1: K = 1
            if P < 1: P = 1
    else:
        K = None
        P = None
    K, P = cOmm.bcast([K, P], root=mAster_rank)

    mesh = random_3D_mesh_of_elements_around(K, **kwargs)

    space = random_3D_space_of_overall_degrees_around(P)

    return mesh, space




def random_3D_FormCaller_of_total_load_around(*args, **kwargs):
    """A wrapper of `random_3D_mesh_and_space_of_total_load_around` and we use the outputs to make a
    3D FormCaller instance."""
    mesh, space = random_3D_mesh_and_space_of_total_load_around(*args, **kwargs)
    return FormCaller(mesh, space)





if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\TESTS\random_objects.py
    # random_3D_mesh_of_elements_around(1)
    random_3D_mesh_and_space_of_total_load_around(100)
