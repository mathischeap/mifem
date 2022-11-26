# -*- coding: utf-8 -*-



from root.config.main import RANK, MASTER_RANK, COMM
from components.freeze.main import FrozenOnly

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np




class _3dCSCG_2Trace_Visualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return self.matplot(**kwargs)

    def matplot(self, density=10000, i=None,
                         colormap='RdBu',
                         num_color_bar_ticks=5):
        """

        :param density:
        :param i: Plot which trace elements?
        :param colormap:
        :param num_color_bar_ticks:
        :return:
        """
        mesh = self._tf_.mesh
        density = int(np.sqrt(density/mesh.trace.elements.GLOBAL_num)) + 1
        xi = eta = sigma = np.linspace(-1, 1, density)
        xyz, v = self._tf_.reconstruct(xi, eta, sigma, i=i)
        xyz = COMM.gather(xyz, root=MASTER_RANK)
        v = COMM.gather(v, root=MASTER_RANK)
        if RANK != MASTER_RANK: return

        XYZ = list()
        V = list()
        for _xyz_, _v_ in zip(xyz, v):
            for i in _xyz_:
                xyz_i = _xyz_[i]
                v_i = _v_[i][0]

                XYZ.append(xyz_i)
                V.append(v_i)

        V = np.array(V)
        del xyz, v

        MAX = np.max(V)
        MIN = np.min(V)
        if MAX == MIN:
            MAX += 0.0001

        bounds = MAX - MIN
        V = V - MIN
        V = V / bounds

        ticks = np.linspace(MAX, MIN, num_color_bar_ticks)

        cmap = getattr(cm, colormap)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(XYZ):
            x, y, z = xyz
            v = V[i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks))
        cb = plt.colorbar(mappable, ax=ax, # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,# extend='min',
                          orientation='vertical', )
        # cb.set_label(r'$ \log_{10}\left( \left| \lambda^h-\varphi_{\mathrm{exact}} \right| \right) $',
        #     labelpad=10, size=15)
        cb.ax.tick_params()#labelsize=13.5)

        plt.show()