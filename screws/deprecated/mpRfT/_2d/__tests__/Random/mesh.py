# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import RANK, MASTER_RANK, COMM
import random
from tests.objects.CSCG._2d.randObj.mesh import random_mesh_of_elements_around as cscg_random_mesh2
from objects.mpRfT._2d.mesh.main import mpRfT2_Mesh

def rm(elements_num,
       N_range = (1,3),
       refinement_intensity=0.5,
       refinement_levels=3,
       exclude_periodic=False,
       domain_boundary_distribution_regularities=None,
       mesh_boundary_num = None,
       mesh_pool=None,
       EDM_pool = None,
    ):
    """

    Parameters
    ----------
    elements_num : int
        The mesh will be of around this many level-0-cells (cscg mesh elements).
    N_range :
    refinement_intensity : float, optional
    refinement_levels : int, optional
    exclude_periodic : bool, optional
    domain_boundary_distribution_regularities
    mesh_boundary_num
    mesh_pool
    EDM_pool

    Returns
    -------

    """
    if RANK == MASTER_RANK:
        dN = random.randint(*N_range)
    else:
        dN = None
    dN = COMM.bcast(dN, root=MASTER_RANK)
    elements_num = int(elements_num / dN) + 10

    cscg = cscg_random_mesh2(elements_num,
           exclude_periodic=exclude_periodic,
           domain_boundary_distribution_regularities=domain_boundary_distribution_regularities,
           mesh_boundary_num = mesh_boundary_num,
           mesh_pool=mesh_pool,
           EDM_pool = EDM_pool)

    #------ Now we randomly generate some refinements -----------------------------------------
    rfd = dict()

    num_elements = cscg.elements.num # local num of mesh-elements (level-0-cells)
    num = int(num_elements * refinement_intensity)

    refine_pool = [str(_) for _ in cscg.elements.indices]
    _2b_refined = random.sample(refine_pool, num)

    lv = 0
    while num > 0 and lv < refinement_levels:

        _2b_p_refined = random.sample(_2b_refined, int(0.1 * random.randint(3,6) * len(_2b_refined)))

        _2b_h_refined = list()
        for c in _2b_refined:
            if c in _2b_p_refined:
                rfd[c] = random.randint(*N_range)
            else:
                _2b_h_refined.append(c)

        _2b_refined = list()
        for c in _2b_h_refined:
            if '-' not in c:
                _2b_refined.extend([c+'-'+str(_) for _ in range(4)])
            else:
                _2b_refined.extend([c+str(_) for _ in range(4)])

        lv += 1

        num = int(len(_2b_refined) * refinement_intensity)

    for c in _2b_refined:
        rfd[c] = random.randint(*N_range)

    mesh = mpRfT2_Mesh(cscg, dN, rfd)

    for i in mesh:
        cell = mesh[i]
        assert cell.N is not None # checking all root cells are made through cells_dict

    return mesh






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/__tests__/Random/mesh.py
    mesh = rm(50)
    mesh.visualization()

