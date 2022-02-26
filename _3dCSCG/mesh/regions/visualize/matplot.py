

from root.config import *
from screws.frozen import FrozenOnly
import matplotlib.pyplot as plt
from matplotlib import cm


class _3dCSCG_Regions_Visualize_Matplot_(FrozenOnly):
    def __init__(self, visualize):
        self._regions_ = visualize._regions_
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.topology(*args, **kwargs)


    def topology(self, aspect='equal'):
        """"""
        if rAnk != mAster_rank: return # only need the master core.

        region_coordinates_pool, connection_pool = self._regions_.___PRIVATE_parse_topology_1___()


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

        H = 0.5
        x_lim, y_lim, z_lim = [list() for _ in range(3)]

        MAP = self._regions_.map
        BNS = self._regions_._domain_.boundaries.names
        colors = cm.get_cmap('cool_r', len(BNS))
        COLORS = dict()
        for i, _ in enumerate(BNS):
            COLORS[_] = colors(i)

        for rn in region_coordinates_pool:
            ax.scatter(*region_coordinates_pool[rn], marker='s')

            ax.text(*region_coordinates_pool[rn], rn,
                    ha='center', va='center', ma='center')

            x, y, z = region_coordinates_pool[rn]

            line = ([x-H, x-H, x-H, x-H, x-H],
                    [y-H, y+H, y+H, y-H, y-H],
                    [z-H, z-H, z+H, z+H, z-H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            line = ([x+H, x+H, x+H, x+H, x+H],
                    [y-H, y+H, y+H, y-H, y-H],
                    [z-H, z-H, z+H, z+H, z-H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            line = ([x-H, x+H],
                    [y-H, y-H],
                    [z-H, z-H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            line = ([x-H, x+H],
                    [y+H, y+H],
                    [z-H, z-H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            line = ([x-H, x+H],
                    [y-H, y-H],
                    [z+H, z+H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            line = ([x-H, x+H],
                    [y+H, y+H],
                    [z+H, z+H])
            ax.plot(*line, '-', color=(0.5,0.5,0.5,0.15))

            SIDES = MAP[rn]
            for s, what in zip('NSWEBF', SIDES):
                if what in BNS:
                    if s == 'N':
                        ax.text(x-H, y, z, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    elif s == 'S':
                        ax.text(x+H, y, z, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    elif s == 'W':
                        ax.text(x, y-H, z, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    elif s == 'E':
                        ax.text(x, y+H, z, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    elif s == 'B':
                        ax.text(x, y, z-H, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    elif s == 'F':
                        ax.text(x, y, z+H, '<'+what+'>', color=COLORS[what], ha='center', va='center', ma='center')
                    else:
                        raise Exception()
                else:
                    assert what[:2] == 'R:'

            if aspect == 'equal':
                x_lim.extend([x-H, x+H])
                y_lim.extend([y-H, y+H])
                z_lim.extend([z-H, z+H])

        for cn in connection_pool:
            ax.plot(*connection_pool[cn], '--', color='gray')

        if aspect == 'equal':
            ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)))

        ax.tick_params(labelsize=10)
        ax.set_xlabel(r'$x$', fontsize=12)
        ax.set_ylabel(r'$y$', fontsize=12)
        ax.set_zlabel(r'$z$', fontsize=12)
        plt.title(self._regions_._domain_.name +
                  ', ID: ' +
                  self._regions_._domain_.parameters['ID'] +
                  ', <regions topology>')
        fig.tight_layout()
        plt.show()
        return fig