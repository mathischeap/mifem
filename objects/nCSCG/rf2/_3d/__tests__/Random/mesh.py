# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 6:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import random
from objects.CSCG._3d.__tests__.Random.mesh import random_mesh_of_elements_around as cscg_random_mesh3
from objects.nCSCG.rf2._3d.mesh.main import _3nCSCG_RF2_Mesh

def random_mesh_of_elements_around(elements_num,
                                   refinement_intensity=0.25,
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
    cscg = cscg_random_mesh3(elements_num,
           exclude_periodic=exclude_periodic,
           domain_boundary_distribution_regularities=domain_boundary_distribution_regularities,
           mesh_boundary_num = mesh_boundary_num,
           mesh_pool=mesh_pool,
           EDM_pool = EDM_pool)

    mesh = _3nCSCG_RF2_Mesh(cscg)

    #------ Now we randomly generate some refinements -----------------------------------------
    num_elements = cscg.elements.num # local num of mesh-elements (level-0-cells)
    refinement_level_num = dict() # how many cells to be refined on each level.
    for l in range(refinement_levels):
        num_elements = int(num_elements * refinement_intensity)
        refinement_level_num[l] = num_elements

    refine_pool = cscg.elements.indices
    for l in refinement_level_num:
        if refinement_level_num[l] == 0: break # no need to do further refinement.

        cells_2b_refined = random.sample(refine_pool, refinement_level_num[l])
        for i in cells_2b_refined:
            mesh.do.refine(i)

        refine_pool = list()
        for i in cells_2b_refined:
            cell = mesh(i)
            for j in range(8): # make this range to be 4 for 2d nCSCG_RF2 mesh
                sub_cell = cell.sub_cells[j]
                refine_pool.append(sub_cell.indices)

    return mesh

if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_3d/__tests__/Random/mesh.py
    mesh = random_mesh_of_elements_around(100)
    mesh.visualize()
