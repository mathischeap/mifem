# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 6:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, cOmm, mAster_rank, np


class _2nCSCG_MeshIDS_Vector_ReGionW(FrozenOnly):
    """data are arranged Region-Wise."""

    def __init__(self, vector):
        """"""
        # below is an over-powerful constraint, we only need region-wise homogeneous, but below is globally homogeneous.
        assert vector.distribution == 'homogeneous', \
            f"distribution={vector.distribution} wrong, only accept homogeneous distribution."
        assert vector.ndim == 2, f"Only 2-dim data can be BCW arranged."
        self._vector_ = vector
        self.___Pr_stack_to_regions___(vector.BCW._BCW_) # signature check will be done in BCW.
        self._freeze_self_()


    def ___Pr_stack_to_regions___(self, BCW):
        """"""
        BCW = cOmm.gather(BCW, root=mAster_rank)
        mesh = self._vector_.mesh
        if rAnk == mAster_rank:
            ___ = dict()
            for _ in BCW:
                ___.update(_)
            BCW = ___

            for _ in BCW: assert np.ndim(BCW[_]) == 3
            assert len(BCW) == mesh.cscg._num_total_elements_

            _sd_ = dict()
            ij = np.shape(BCW[0][0])
            I, J = ij
            ALL_element_global_numbering_ = \
                mesh.cscg.___PRIVATE_generate_ALL_element_global_numbering___()

            for Rn in ALL_element_global_numbering_:
                region_data_shape = [ij[i] * mesh.cscg._element_layout_[Rn][i] for i in range(2)]
                _sd_[Rn] = [np.zeros(region_data_shape), np.zeros(region_data_shape)]
                for j in range(mesh.cscg._element_layout_[Rn][1]):
                    for i in range(mesh.cscg._element_layout_[Rn][0]):
                        _sd_[Rn][0][i * I:(i + 1) * I, j * J:(j + 1) * J] = \
                            BCW[ALL_element_global_numbering_[Rn][i, j]][0]
                        _sd_[Rn][1][i * I:(i + 1) * I, j * J:(j + 1) * J] = \
                            BCW[ALL_element_global_numbering_[Rn][i, j]][1]

            self._RGW_ = _sd_

        else:
            self._RGW_ = None

    def __getitem__(self, rn):
        assert rAnk == mAster_rank, f"please access RGW only through master core."
        return self._RGW_[rn]

    def __iter__(self):
        assert rAnk == mAster_rank, f"please access RGW only through master core."
        for rn in self._RGW_:
            yield rn



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
