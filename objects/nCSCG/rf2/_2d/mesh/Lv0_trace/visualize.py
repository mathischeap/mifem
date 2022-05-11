# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 10:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import np, rAnk, mAster_rank, cOmm

import matplotlib.pyplot as plt
import matplotlib


class _2nCSCG_Lv0TraceVisualize(FrozenOnly):
    """"""

    def __init__(self, lv0trace):
        """"""
        self._trace_ = lv0trace
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)


    def matplot(self, density=5, show_segments=False,
        saveto=None, usetex=False):
        """"""
        xi = np.linspace(-1, 1, density)
        all_trace_segments = self._trace_.elements.segments

        SEG_DATA = dict()

        for i in all_trace_segments: # go through all local lv0-trace-elements (cscg trace elements)
            segments = all_trace_segments[i]
            for segment in segments:
                xy = segment.coordinate_transformation.mapping(xi)
                SEG_DATA[segment.__repr__()] = xy

        SEG_DATA = cOmm.gather(SEG_DATA, root=mAster_rank)

        if rAnk != mAster_rank: return

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
        plt.close()





if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/Lv0_trace/visualize.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)

    # from root.read.main import read
    # mesh = read('test_mesh.mi')

    trace = mesh.Lv0trace
    trace.elements.do.update()
    # print(trace.elements.segments)
    trace.visualize()
    mesh.visualize()