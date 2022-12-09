# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly

import matplotlib.pyplot as plt
from matplotlib import cm
from root.config.main import *

class _3dCSCG_ScalarField_matplot_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.boundary_values(*args, **kwargs)

    def boundary_values(self,
        density=5000, colormap='coolwarm',
        num_color_bar_ticks=5):
        """"""
        mesh = self._f_.mesh
        NUM = mesh.trace.elements.global_num
        density = int((density/NUM)**0.5 + 1)

        x = y = z = np.linspace(-1, 1, density)

        xyz, v = self._f_.reconstruct(x, y, z, where='trace-element')

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
        cb.ax.tick_params(labelsize=13.5)


        plt.show()
