


import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly

from root.config.main import rAnk, mAster_rank, np
import matplotlib.pyplot as plt



class _3dCSCG_Mesh_Boundary_VIS(FrozenOnly):
    """"""
    def __init__(self, boundary):
        """"""
        self._boundary_ = boundary
        self._mesh_ = boundary.mesh
        self._default_ = 'matplot'
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, density=1000, aspect='equal',):
        """"""
        region_sides = self._boundary_.region_sides

        if rAnk != mAster_rank: return

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

        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        for rn in domain.regions:
            region = domain.regions[rn]
            for i, sn in enumerate('NSWEBF'):
                is_boundary = regions.sides_on_domain_boundaries[rn][i]
                rs = region.sides[sn]
                xyz = rs.coordinate_transformation.mapping(rrr, sss)
                if aspect == 'equal':
                    x, y, z = xyz
                    x_lim.append(np.max(x))
                    x_lim.append(np.min(x))
                    y_lim.append(np.max(y))
                    y_lim.append(np.min(y))
                    z_lim.append(np.max(z))
                    z_lim.append(np.min(z))
                if is_boundary:
                    ax.plot_surface(*xyz, color=(1,1,1,0.25))

        ONES = np.ones(density)
        SPACING = np.linspace(0,1, density)
        for region_side in region_sides:
            # print(bn, region_side)
            rn, side = region_side.split('-')
            region = domain.regions[rn]
            rs = region.sides[side]
            xyz = rs.coordinate_transformation.mapping(rrr, sss)
            ax.plot_surface(*xyz, color=(1,0.95,1,0.5))

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
            for i in range(Lr+1):
                lines.append([r[i]*ONES, SPACING])
            for i in range(Ls+1):
                lines.append([SPACING, s[i]*ONES])

            for line in lines:
                xyz = rs.coordinate_transformation.mapping(*line)
                ax.plot(*xyz, color='k', linewidth=0.8)

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        plt.title(domain.name + ', ID: '+
                  domain.parameters['ID'] +
                  f', <mesh boundary: {self._boundary_.name}>')

        fig.tight_layout()
        plt.show()
        plt.close(fig)





if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\mesh\boundaries\boundary\visualize.py
    from objects.CSCG._3d.master import MeshGenerator
    mesh = MeshGenerator('bridge_arch_cracked')([2,3,5])
    # mesh = MeshGenerator('crazy_periodic')([4,2,[1,2,4,2,1]])
    boundaries = mesh.boundaries

    # boundary = boundaries['Bottom']
    # boundary = boundaries['Crack']
    boundary = boundaries['Left_Wall']
    # boundary = boundaries['Right_Floor']
    # boundary = boundaries['Front']

    boundary.visualize()