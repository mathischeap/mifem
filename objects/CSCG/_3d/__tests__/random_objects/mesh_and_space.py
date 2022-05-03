# -*- coding: utf-8 -*-
from root.config.main import rAnk, mAster_rank, cOmm
import random
from objects.CSCG._3d.__tests__.random_objects.mesh import random_mesh_of_elements_around
from objects.CSCG._3d.__tests__.random_objects.space import random_space_of_degrees_around


def random_mesh_and_space_of_total_load_around(total_load, **kwargs):
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

    mesh = random_mesh_of_elements_around(K, **kwargs)

    space = random_space_of_degrees_around(P)

    return mesh, space