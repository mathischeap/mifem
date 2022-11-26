# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK, np

import matplotlib.pyplot as plt
from matplotlib import cm


class _3dCSCG_MeshBoundaries_VIS(FrozenOnly):
    """"""
    def __init__(self, bds):
        """"""
        self._boundaries_ = bds
        self._mesh_ = bds._mesh_
        self._freeze_self_()



    def __call__(self, density=1000, aspect='equal',):
        """"""
        # we can do everything in the master core.

        boundaries_name = self._boundaries_.names
        boundaries_num = len(boundaries_name)
        range_of_region_sides = self._boundaries_.range_of_region_sides

        if RANK != MASTER_RANK: return

        domain = self._mesh_.domain

        density = int(np.ceil(np.sqrt(density/(domain.regions.num*6))))
        rrr = sss = np.linspace(0,1,density)
        rrr, sss = np.meshgrid(rrr, sss, indexing='ij')

        regions = domain.regions
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

        if boundaries_num > 10:
            corlormap = 'viridis'
        else:
            corlormap = 'tab10'

        color = cm.get_cmap(corlormap, boundaries_num)
        colors = []
        boundary_name_color_dict = dict()
        for j in range(boundaries_num):
            temp_color = list(color(j))
            temp_color[-1] = 0.5
            colors.append(temp_color)
        for j, bn in enumerate(boundaries_name):
            boundary_name_color_dict[bn] = colors[j]

        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for rn in domain.regions:
            region = domain.regions[rn]
            for i, sn in enumerate('NSWEBF'):
                is_boundary = regions.sides_on_domain_boundaries[rn][i]

                rs = region.sides[sn]

                xyz = rs.coordinate_transformation.mapping(rrr, sss)

                # find the ratio.
                if aspect == 'equal':
                    x, y, z = xyz
                    x_lim.append(np.max(x))
                    x_lim.append(np.min(x))
                    y_lim.append(np.max(y))
                    y_lim.append(np.min(y))
                    z_lim.append(np.max(z))
                    z_lim.append(np.min(z))

                if is_boundary:
                    ax.plot_surface(*xyz, color=(1,1,1,0.4))

        ONES = np.ones(density)
        SPACING = np.linspace(0,1, density)
        for bn in range_of_region_sides:

            for region_side in range_of_region_sides[bn]:
                # print(bn, region_side)
                rn, side = region_side.split('-')
                region = domain.regions[rn]
                rs = region.sides[side]
                xyz = rs.coordinate_transformation.mapping(rrr, sss)
                ax.plot_surface(*xyz, color=boundary_name_color_dict[bn])

                x, y, z = xyz
                x_range = [np.min(x), np.max(x)]
                y_range = [np.min(y), np.max(y)]
                z_range = [np.min(z), np.max(z)]
                x_mid, y_mid, z_mid = np.mean(x_range), np.mean(y_range), np.mean(z_range)

                ax.text(x_mid, y_mid, z_mid, bn,
                        color='k',
                        ha='center', va='center', ma='center')

                spacing = self._mesh_.elements.spacing[rn]
                layout = self._mesh_.elements.layout[rn]
                if side in 'NS':
                    r, s = spacing[1], spacing[2]
                    Lr, Ls = layout[1], layout[2]
                elif side in 'WE':
                    r, s = spacing[0], spacing[2]
                    Lr, Ls = layout[0], layout[2]
                elif side in 'BF':
                    r, s = spacing[0], spacing[1]
                    Lr, Ls = layout[0], layout[1]
                else:
                    raise Exception()

                lines = list()
                for i in range(Lr + 1):
                    lines.append([r[i] * ONES, SPACING])
                for i in range(Ls + 1):
                    lines.append([SPACING, s[i] * ONES])

                for line in lines:
                    xyz = rs.coordinate_transformation.mapping(*line)
                    ax.plot(*xyz, color='k', linewidth=1)

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        plt.title(domain.name + ', ID: '+
                  domain.parameters['ID'] +
                  ', <domain>')

        fig.tight_layout()
        plt.show()
        plt.close(fig)





if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\mesh\boundaries\visualize.py
    from objects.CSCG._3d.master import MeshGenerator
    mesh = MeshGenerator('bridge_arch_cracked')([4,2,[1,2,4,2,1]])
    # mesh = MeshGenerator('crazy_periodic')([4,2,[1,2,4,2,1]])
    boundaries = mesh.boundaries

    # mesh.domain.visualize(
    #     show_internal_region_sides=True,
    #     show_boundary_names=False,
    #     distinguish_boundaries_by_color=False,)

    boundaries.visualize()