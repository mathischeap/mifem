

import sys
if './' not in sys.path: sys.path.append('visualize/')
from screws.freeze.base import FrozenOnly

from root.config.main import rAnk, mAster_rank, np, sEcretary_rank
import matplotlib.pyplot as plt

class _3dCSCG_MeshElements_VIS(FrozenOnly):
    """"""
    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()

    def __call__(self, density=1000, aspect='equal',):
        """We plot the local mesh elements in the whole computational domain.

        :return:
        """

        # if rAnk != mAster_rank: return
        domain = self._elements_._mesh_.domain
        density = int(np.ceil(np.sqrt(density/(domain.regions.num*6))))
        if density <= 2: density = 3

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

        #-------- plot the computational domain --------------------------------
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
                    ax.plot_surface(*xyz, color=(1,1,1,0.3))

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        # ----- plot all local mesh elements -------------------------------------
        ELE = self._elements_.indices
        r = np.linspace(-1,1, density)
        O = np.ones(density)
        WB = (r, -1*O, -1*O)
        EB = (r, +1*O, -1*O)
        WF = (r, -1*O, +1*O)
        EF = (r, +1*O, +1*O)
        NB = (-1*O, r, -1*O)
        SB = (+1*O, r, -1*O)
        NF = (-1*O, r, +1*O)
        SF = (+1*O, r, +1*O)
        NW = (-1*O, -1*O, r)
        NE = (-1*O, +1*O, r)
        SW = (+1*O, -1*O, r)
        SE = (+1*O, +1*O, r)

        LINES = dict()
        for i in ELE:
            LINES[i] = dict()
            element = self._elements_[i]
            LINES[i]['WB'] = element.coordinate_transformation.mapping(*WB)
            LINES[i]['EB'] = element.coordinate_transformation.mapping(*EB)
            LINES[i]['WF'] = element.coordinate_transformation.mapping(*WF)
            LINES[i]['EF'] = element.coordinate_transformation.mapping(*EF)
            LINES[i]['NB'] = element.coordinate_transformation.mapping(*NB)
            LINES[i]['SB'] = element.coordinate_transformation.mapping(*SB)
            LINES[i]['NF'] = element.coordinate_transformation.mapping(*NF)
            LINES[i]['SF'] = element.coordinate_transformation.mapping(*SF)
            LINES[i]['NW'] = element.coordinate_transformation.mapping(*NW)
            LINES[i]['NE'] = element.coordinate_transformation.mapping(*NE)
            LINES[i]['SW'] = element.coordinate_transformation.mapping(*SW)
            LINES[i]['SE'] = element.coordinate_transformation.mapping(*SE)

        for i in LINES:
            for edge in LINES[i]:
                plt.plot(*LINES[i][edge], 'k', linewidth=0.6)

        # ---------------------------------------------------------------------------------
        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=15)
        ax.set_ylabel(r'$y$', fontsize=15)
        ax.set_zlabel(r'$z$', fontsize=15)
        if rAnk == mAster_rank:
            plt.title(domain.name + f'; mesh elements in the MASTER core #{rAnk}')
        elif rAnk == sEcretary_rank:
            plt.title(domain.name + f'; mesh elements in the SECRETARY core #{rAnk}')
        else:
            plt.title(domain.name + f'; mesh elements in core #{rAnk}')

        fig.tight_layout()
        plt.show()
        plt.close(fig)




if __name__ == '__main__':
    # mpiexec -n 3 python _3dCSCG\mesh\elements\visualize\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 3, 3]
    # mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)


    mesh.elements.visualize()