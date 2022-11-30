# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 3:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from root.config.main import RANK, MASTER_RANK, COMM

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


class _3dCSCG_2LocalTrace_Visualize_Matplot(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.boundary_value(*args, **kwargs)

    def boundary_value(self,
                       density=10000,
                       colormap='RdBu',
                       num_color_bar_ticks=5):
        """

        Parameters
        ----------
        density
        colormap
        num_color_bar_ticks

        Returns
        -------

        """
        mesh = self._ltf_.mesh
        density = int(np.sqrt(density/mesh.trace.elements.GLOBAL_num)) + 1
        xi = eta = sigma = np.linspace(-1, 1, density)

        range_element_side = mesh.boundaries.range_of_element_sides

        local_elements = list()
        local_sides = list()
        for bn in range_element_side:
            elements_sides = range_element_side[bn]
            for element_side in elements_sides:
                element = int(element_side[:-1])
                side = element_side[-1]
                local_elements.append(element)
                local_sides.append(side)

        local_element_set = set(local_elements)

        xyz, v = self._ltf_.reconstruct(xi, eta, sigma, element_range=local_element_set, ravel=False)

        xyz = COMM.gather(xyz, root=MASTER_RANK)
        v = COMM.gather(v, root=MASTER_RANK)
        local_elements = COMM.gather(local_elements, root=MASTER_RANK)
        local_sides = COMM.gather(local_sides, root=MASTER_RANK)
        if RANK != MASTER_RANK: return

        ___ = dict()
        for _ in xyz:
            ___.update(_)
        xyz = ___

        ___ = dict()
        for _ in v:
            ___.update(_)
        v = ___

        plot_xyz = list()
        plot_v = list()
        for elements, sides in zip(local_elements, local_sides):
            for element, side in zip(elements, sides):
                plot_xyz.append(xyz[element][side])
                plot_v.append(v[element][side][0])

        plot_v = np.array(plot_v)
        del xyz, v

        MAX = np.max(plot_v)
        MIN = np.min(plot_v)
        if MAX == MIN:
            MAX += 0.0001

        bounds = MAX - MIN
        plot_v = plot_v - MIN
        plot_v = plot_v / bounds

        ticks = np.linspace(MAX, MIN, num_color_bar_ticks)

        cmap = getattr(cm, colormap)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(plot_xyz):
            x, y, z = xyz
            v = plot_v[i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks))
        cb = plt.colorbar(
            mappable,
            ax=ax, # ticks=np.linspace(0,1,num_ticks),
            shrink=1,
            aspect=20, # extend='min',
            orientation='vertical',
        )
        cb.set_label(
            f'{self._ltf_.name}',
            labelpad=10, size=15
        )
        cb.ax.tick_params()#labelsize=13.5)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        plt.show()




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
