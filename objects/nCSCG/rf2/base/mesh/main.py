# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 2:10 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import rAnk, mAster_rank, cOmm
from screws.freeze.main import FrozenClass
from objects.nCSCG.rf2.base.mesh.historic.main import nCSCG_Historic



class nCSCG_RF2_MeshBase(FrozenClass):
    """"""
    def __init__(self, cscg):
        """We build a non-conforming mesh based on a cscg mesh.

        So this cscg mesh will be the coarsest case of this non-conforming mesh.

        Parameters
        ----------
        cscg
        """
        assert cscg.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'), f"I need a CSCG mesh."
        self._cscg_ = cscg
        self.___base_mesh_cells___ = None
        self._space_ = None # will be set when we set the mesh of the space
        self._historic_ = None
        self._locker_ = True
        self._signature_ = None # will be updated when updating mesh

    @property
    def signature(self):
        """will be updated when updating mesh"""
        return self._signature_

    def __call__(self, indices, dynamic=False):
        """

        Parameters
        ----------
        indices
        dynamic : bool, optional
            If it is True, we will refine it when the cell is not refined yet. Otherwise, we will
            raise exception.

        Returns
        -------

        """
        if not isinstance(indices, tuple):
            indices = (indices,)
        else:
            pass

        cells = self.___base_mesh_cells___
        for _, i in enumerate(indices):
            if cells is None:
                if dynamic: # refine it when the cell is not refined yet; not safe; for test purposes.
                    cell.do.refine()
                    cells = cell.sub_cells
                else:
                    raise Exception(f"cell[{indices[:_]}] has no level: >{_}< sub-cells.")
            else:
                pass
            cell = cells[i]
            cells = cell.sub_cells

        return cell

    def __iter__(self):
        """iter through all indices of local smallest cells."""
        cells = self.___base_mesh_cells___
        for i in cells:
            cell = cells[i]
            for j in cell:
                yield j

    def __contains__(self, item):
        raise Exception(f"__contains__ of nCSCG mesh is disabled.")

    @property
    def ___parameters___(self):
        """
        This `parameters` is used to compare if meshes are the same. Therefore, the
        `___parameters___` should uniquely identify a mesh. We also use it tor save and restore a mesh.

        So it is mandatory for saving a mesh.
        """
        cscg_parameters = self._cscg_.standard_properties.parameters
        parameters = dict()
        parameters['type'] = self.__class__.__name__
        parameters['cscg'] = cscg_parameters
        parameters['space'] = self.space._PRM

        #Note that the history of the refinement will be lost!
        refinements = dict()
        for ind in self:
            if len(ind) > 1: # this level-0-cell is refined.
                mesh_element = ind[0]
                region_indices = self.cscg.do.find.region_name_and_local_indices_of_element(mesh_element)
                region, indices = region_indices
                key = region + '|' + str(indices)[1:-1]
                if key not in refinements:
                    refinements[key] = list()
                refinements[key].append(ind[1:])
            else:
                pass

        refinements = cOmm.gather(refinements, root=mAster_rank)
        if rAnk == mAster_rank:
            _R_ = dict()
            for _r in refinements:
                _R_.update(_r)
            refinements = _R_
        parameters['refinements'] = refinements # All info in master, None in slaves.

        return parameters

    def __eq__(self, other):
        """Regardless of the history of a mesh, we only compare they current conditions."""
        return self.standard_properties.parameters == other.standard_properties.parameters

    @property
    def cscg(self):
        return self._cscg_

    @property
    def ndim(self):
        return self._cscg_.ndim

    @property
    def domain(self):
        """"""
        return self._cscg_.domain

    @property
    def base_cells(self):
        """The base cells, namely, the level-0-cells, namely, the cscg mesh elements.

        These cells cannot be deleted.
        """
        return self.___base_mesh_cells___

    @property
    def historic(self):
        """"""
        if self._historic_ is None:
            self._historic_ = nCSCG_Historic(self)
        return self._historic_

    @property
    def space(self):
        return self._space_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/base/mesh/main.py
    from  objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    # # from  objects.nCSCG.rf2._3d.__tests__.Random.mesh import random_mesh_of_elements_around as rm3
    # #
    mesh2 = rm2(100)
    # # mesh3 = rm3(200)
    # #
    mesh2.visualize()
    # # print(mesh2.standard_properties.parameters)
    # #
    # from root.save import save
    # #
    # save(mesh2, '___test_mesh2_sr__')

    # from root.read.main import read
    # #
    # mesh = read('___test_mesh2_sr__')
    # mesh.visualize()