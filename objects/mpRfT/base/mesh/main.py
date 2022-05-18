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


class mpRfT_MeshBase(FrozenClass):
    """"""
    def __init__(self, cscg, dN, cells_dict):
        """

        Parameters
        ----------
        cscg
        dN
        cells_dict
        """
        assert cscg.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'), f"I need a CSCG mesh."
        self.___cscg___ = cscg
        assert isinstance(dN, int) and dN > 0, f"dN={dN} wrong, must be a positive integer."
        self.___dN___ = dN
        self.___basic_cells___ = None
        self.___Pr_initialize_cells___(dN, cells_dict)

    def ___Pr_initialize_cells___(self, dN, cells_dict):
        """"""
        raise NotImplementedError()

    @property
    def dN(self):
        return self.___dN___

    @property
    def basic_cells(self):
        """The basic cells, namely, the level-0-cells, namely, the cscg mesh elements.

        These cells cannot be diluted.
        """
        return self.___basic_cells___

    def __getitem__(self, indices):
        """

        Parameters
        ----------
        indices

        Returns
        -------

        """
        if isinstance(indices, int):
            indices = (indices,)
        else:
            assert isinstance(indices, tuple), f"put indices of cell in a tuple."

        cells = self.basic_cells
        for _, i in enumerate(indices):
            if cells is None:
                raise Exception(f"cell[{indices[:_]}] has no level: >{_}< sub-cells.")
            else:
                pass
            cell = cells[i]
            cells = cell.sub_cells

        return cell

    def __iter__(self):
        """iter through all indices of local root cells."""
        cells = self.basic_cells
        for i in cells:
            cell = cells[i]
            for j in cell:
                yield j

    def __contains__(self, item):
        raise Exception(f"`__contains__` method of mpRfT mesh is disabled.")

    def __len__(self):
        raise Exception(f"`__len__` method of mpRfT mesh is disabled.")

    def __call__(self, item):
        raise Exception(f"`__call__` method of mpRfT mesh is disabled.")

    @property
    def ___parameters___(self):
        """It is mandatory for saving a mesh.
        """
        cscg_parameters = self.___cscg___.standard_properties.parameters
        parameters = dict()
        parameters['type'] = self.__class__.__name__
        parameters['cscg'] = cscg_parameters
        parameters['dN'] = self.dN

        #Note that the history of the refinement will be lost!
        refinements = dict()
        for ind in self:
            cell = self[ind]

            if len(ind) > 1:
                mesh_element = ind[0]
                region_indices = self.cscg.do.find.region_name_and_local_indices_of_element(mesh_element)
                region, indices = region_indices
                key = region + '|' + str(indices)[1:-1]
                if key not in refinements:
                    refinements[key] = dict()

                rp = cell.__repr__()
                assert cell.N is not None, f"cell:[{rp}] N is None."
                rp = rp.split('-')[1]
                refinements[key][rp] = cell.N

            else:
                cell = self[ind]
                if cell.N != self.dN: # p-refined basic-root-cell, store it.


                    mesh_element = ind[0]
                    region_indices = self.cscg.do.find.region_name_and_local_indices_of_element(mesh_element)
                    region, indices = region_indices
                    key = region + '|' + str(indices)[1:-1]
                    if key not in refinements:
                        refinements[key] = dict()

                    assert cell.N is not None, f"cell:[{cell}] N is None."
                    assert len(refinements[key]) == 0, f"must be the case."
                    refinements[key][''] = cell.N

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
        return self.___cscg___

    @property
    def ndim(self):
        return self.___cscg___.ndim







if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/base/mesh/main.py
    pass