# -*- coding: utf-8 -*-
from objects.CSCG._3d.master import SpaceInvoker
from root.config.main import RANK, MASTER_RANK, np, COMM
import random


def random_space_of_degrees_around(overall_degrees):
    """

    :param overall_degrees: p0 * p1 * p2
    :return:
    """
    assert isinstance(overall_degrees, (int, float)) and overall_degrees >= 1, \
        f"overall_degrees = {overall_degrees} is wrong, must be a number and >= 1."
    if RANK == MASTER_RANK:

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
    space_name, DIS = COMM.bcast([space_name, DIS], root=MASTER_RANK)

    space = SpaceInvoker(space_name)(DIS)

    return space