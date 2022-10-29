# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 4:25 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import RANK, COMM, MASTER_RANK, np


class mpRfT2_Mesh_rcWds_Scalar_rgW(FrozenOnly):
    """data are arranged Region-Wise."""

    def __init__(self, scalar):
        """"""
        # below is an over-powerful constraint, we only need base-cell-wise uniform, but below is globally uniform.
        assert scalar._distribution_ == 'uniform', \
            f"distribution={scalar._distribution_} wrong, only accept uniform distribution."
        assert scalar._ndim_ == 2, f"Only 2-dim data can be region-wise arranged."
        assert scalar._isfull_, f"only full rcWds scalar can be reshaped into rgW!"
        self._scalar_ = scalar
        self.___Pr_stack_to_regions___(scalar.bcW._BCW_) # signature check will be done in BCW.
        self._freeze_self_()

    def ___Pr_stack_to_regions___(self, BCW):
        """"""
        BCW = COMM.gather(BCW, root=MASTER_RANK)
        mesh = self._scalar_._mesh_
        if RANK == MASTER_RANK:
            ___ = dict()
            for _ in BCW:
                ___.update(_)
            BCW = ___

            for _ in BCW: assert np.ndim(BCW[_]) == 2
            assert len(BCW) == mesh.cscg._num_total_elements_

            _sd_ = dict()
            ij = np.shape(BCW[0])
            I, J = ij
            ALL_element_global_numbering_ = \
                mesh.cscg.___PRIVATE_generate_ALL_element_global_numbering___()

            for Rn in ALL_element_global_numbering_:
                region_data_shape = [ij[i] * mesh.cscg._element_layout_[Rn][i] for i in range(2)]
                _sd_[Rn] = np.zeros(region_data_shape)
                for j in range(mesh.cscg._element_layout_[Rn][1]):
                    for i in range(mesh.cscg._element_layout_[Rn][0]):
                        _sd_[Rn][i * I:(i + 1) * I, j * J:(j + 1) * J] = \
                            BCW[ALL_element_global_numbering_[Rn][i, j]]
            self._RGW_ = _sd_

        else:
            self._RGW_ = None

    def __getitem__(self, rn):
        assert RANK == MASTER_RANK, f"please access RGW only through master core."
        return self._RGW_[rn]

    def __iter__(self):
        assert RANK == MASTER_RANK, f"please access RGW only through master core."
        for rn in self._RGW_:
            yield rn



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
