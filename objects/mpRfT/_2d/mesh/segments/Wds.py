# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/31 10:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import matplotlib.pyplot as plt
from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, mAster_rank, cOmm


class mpRfT2_Mesh_SegmentWiseDataStructure(FrozenOnly):
    """segment-Wise Data Structure."""

    def __init__(self, segments):
        """"""
        self._segments_ = segments
        self._dd_ = None
        self._freeze_self_()


    def __call__(self, data_dict):
        self._dd_ = data_dict
        return self


    def visualization(self, xy):
        """

        Parameters
        ----------
        xy

        Returns
        -------

        """

        xy = cOmm.gather(xy._dd_, root=mAster_rank)
        vv = cOmm.gather(self._dd_, root=mAster_rank)

        if rAnk != mAster_rank: return

        XY = dict()
        for __ in xy:
            XY.update(__)

        VV = dict()
        for __ in vv:
            VV.update(__)


        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection='3d')
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

        for rp in XY:
            ax.plot(*XY[rp], VV[rp], lw=0.75, color='k')

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        plt.show()
        plt.close(fig)








if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
