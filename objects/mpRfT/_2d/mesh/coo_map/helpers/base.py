# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 9:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_CooMapBase(FrozenOnly):
    """"""

    @property
    def ___Pr_is_mpRfT2_mesh_coo_map___(self):
        return True

    @property
    def distribution(self):
        raise NotImplementedError()


    @property
    def ___Pr_rcMC_key___(self):
        """A key implying the value for metric involved computing in each root-cell."""
        raise NotImplementedError()

    def ___Pr_rcMC_nodes___(self, rp):
        raise NotImplementedError()

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
