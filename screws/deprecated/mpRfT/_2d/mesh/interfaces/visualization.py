# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/21 3:10 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import np, RANK, MASTER_RANK, COMM
import matplotlib.pyplot as plt
import matplotlib
from screws.miscellaneous.color_names import tableau_color_names


class mpRfT2_Mesh_Interfaces_Visualization(FrozenOnly):
    """"""

    def __init__(self, IFs):
        """"""
        self._IFs_ = IFs
        self._mesh_ = IFs._mesh_
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    def matplot(self, density=5,
        saveto=None, usetex=False):
        """"""
        mesh = self._mesh_
        xi = np.linspace(-1, 1, density)

        SEG_DATA = dict()

        for sg in self._mesh_.segments: # go through all local segments.
            xy = sg.coordinate_transformation.mapping(xi)
            SEG_DATA[sg.__repr__()] = xy

        IFs = self._IFs_
        IF_keys = list(IFs.___interfaces___.keys())

        SEG_DATA = COMM.gather(SEG_DATA, root=MASTER_RANK)
        IF_keys = COMM.gather(IF_keys, root=MASTER_RANK)

        if RANK != MASTER_RANK: return

        _SD = dict()
        for _s in SEG_DATA:
            _SD.update(_s)
        SEG_DATA = _SD

        _SD = set()
        for _s in IF_keys:
            _SD.update(_s)
        IF_keys = _SD

        if saveto is not None: matplotlib.use('Agg')
        plt.rc('text', usetex=usetex)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        plt.xlabel(r"$x$", fontsize=12)
        plt.ylabel(r"$y$", fontsize=12)
        plt.tick_params(axis='both', which='both', labelsize=12)

        color_names = tableau_color_names()
        color_len = len(color_names)

        for i, IFk in enumerate(IF_keys):
            segments = IFk[3:]
            if '|' in segments:
                segments = segments.split('|')
            else:
                segments = [segments,]

            for sg_rp in segments:
                xy = SEG_DATA.pop(sg_rp)
                x, y = xy
                ax.scatter(x[0], y[0], s=5, color='k', marker='s')
                ax.scatter(x[-1], y[-1], s=2, color='k')
                plt.plot(x, y, linewidth=1.25, color = color_names[i % color_len])

        assert len(SEG_DATA) == 0

        #---------------------- save to --------------------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')

        #=======================================================================================
        plt.close()






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
