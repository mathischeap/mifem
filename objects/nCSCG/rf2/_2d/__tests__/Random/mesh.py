# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import rAnk, mAster_rank, cOmm
import random
from objects.CSCG._2d.__tests__.Random.mesh import random_mesh_of_elements_around as cscg_random_mesh2
from objects.nCSCG.rf2._2d.mesh.main import _2nCSCG_RF2_Mesh

def random_mesh_of_elements_around(elements_num,
       refinement_intensity=0.50,
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
    if rAnk == mAster_rank:
        dN = random.randint(1,3)
    else:
        dN = None
    dN = cOmm.bcast(dN, root=mAster_rank)
    elements_num = int(elements_num / dN) + 10

    cscg = cscg_random_mesh2(elements_num,
           exclude_periodic=exclude_periodic,
           domain_boundary_distribution_regularities=domain_boundary_distribution_regularities,
           mesh_boundary_num = mesh_boundary_num,
           mesh_pool=mesh_pool,
           EDM_pool = EDM_pool)

    mesh = _2nCSCG_RF2_Mesh(cscg, dN)

    #------ Now we randomly generate some refinements -----------------------------------------
    num_elements = cscg.elements.num # local num of mesh-elements (level-0-cells)
    refinement_level_num = dict() # how many cells to be refined on each level.
    for l in range(refinement_levels):
        num_elements = int(num_elements * refinement_intensity)
        refinement_level_num[l] = num_elements

    mesh.do.unlock()
    refine_pool = cscg.elements.indices
    for l in refinement_level_num:
        if refinement_level_num[l] == 0: break # no need to do further refinement.

        cells_2b_refined = random.sample(refine_pool, refinement_level_num[l])
        for i in cells_2b_refined:
            mesh.do.refine(i)

        refine_pool = list()
        for i in cells_2b_refined:
            cell = mesh(i)
            for j in range(4): # make this range to be 8 for 3d nCSCG_RF2 mesh
                sub_cell = cell.sub_cells[j]
                refine_pool.append(sub_cell.indices)

    mesh.do.update()
    return mesh




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/__tests__/Random/mesh.py
    mesh = random_mesh_of_elements_around(300)
    mesh.visualize()
