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
        self._do_ = None
        self._visualize_ = None
        self._freeze_self_()

    def __call__(self, indices, dynamic=False):
        """"""
        if isinstance(indices, int):
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

        refinements = list()
        for ind in self:
            if len(ind) > 1: # this level-0-cell is refined.
                refinements.append(ind)

        refinements = cOmm.gather(refinements, root=mAster_rank)
        if rAnk == mAster_rank:
            _R_ = list()
            for _r in refinements:
                _R_.extend(_r)
            refinements = _R_

        parameters['refinements'] = refinements # All info in master, None in slaves.

        return parameters

    def __eq__(self, other):
        """"""
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





if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/base/mesh/base.py
    pass