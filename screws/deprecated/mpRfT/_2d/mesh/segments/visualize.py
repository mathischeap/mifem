# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/10 10:59 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import np, RANK, MASTER_RANK, COMM

import matplotlib.pyplot as plt
import matplotlib




class mpRfT2_Mesh_Segments_Visualize(FrozenOnly):
    """"""

    def __init__(self, segments):
        """"""
        self._segments_ = segments
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)



    def matplot(self, density=5, show_segments=False,
        saveto=None, usetex=False
        ):
        """"""
        xi = np.linspace(-1, 1, density)

        SEG_DATA = dict()

        for sg in self._segments_: # go through all local lv0-trace-elements (cscg trace elements)
            xy = sg.coordinate_transformation.mapping(xi)
            SEG_DATA[sg.__repr__()] = xy

        SEG_DATA = COMM.gather(SEG_DATA, root=MASTER_RANK)

        if RANK != MASTER_RANK: return

        _SD = dict()
        for _s in SEG_DATA:
            _SD.update(_s)
        SEG_DATA = _SD
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
        for segment in SEG_DATA:
            xy = SEG_DATA[segment]
            x, y= xy
            plt.plot(x, y, linewidth=0.75)
            ax.scatter(x[0], y[0], s=2, color='k')
            ax.scatter(x[-1], y[-1], s=2, color='k')
            if show_segments:
                x = np.mean(x)
                y = np.mean(y)
                plt.text(x, y, segment, ha='center', va='center')

        #---------------------- save to --------------------------------------------------------
        if saveto is None:
            plt.show()
        else:
            plt.savefig(saveto, bbox_inches='tight')

        #=======================================================================================
        plt.close()






if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/segments/visualize.py
    # from objects.nCSCG.rfT2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    # mesh = rm2(100, refinement_intensity=0.5)

    from root.read.main import read
    mesh = read('test_mesh.mi')

    segments = mesh.segments

    segments.visualize(show_segments=True)
