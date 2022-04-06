



from root.config.main import rAnk, mAster_rank, cOmm
from screws.freeze.main import FrozenOnly

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np




class _3dCSCG_Tr_Visualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return self.matplot(**kwargs)

    def matplot(self, density=10000, i=None, colormap='RdBu', levels=5):
        """

        :param density:
        :param i: Plot which trace elements?
        :param colormap:
        :param levels:
        :return:
        """
        mesh = self._Tr_.mesh
        density = int(np.sqrt(density/mesh.trace.elements.GLOBAL_num)) + 1
        xi = eta = sigma = np.linspace(-1, 1, density)
        xyz, v = self._Tr_.reconstruct(xi, eta, sigma, i=i)
        xyz = cOmm.gather(xyz, root=mAster_rank)
        v = cOmm.gather(v, root=mAster_rank)

        if rAnk != mAster_rank: return

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

        ticks = np.linspace(MAX, MIN, levels)

        cmap = getattr(cm, colormap)

        fig = plt.figure(figsize=(12, 8))
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

        plt.close()