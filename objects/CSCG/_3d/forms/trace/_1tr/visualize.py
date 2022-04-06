



from root.config.main import rAnk, mAster_rank, cOmm
from screws.freeze.main import FrozenOnly

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np




class _3dCSCG_1Trace_Visualize(FrozenOnly):
    """The visualization property/component of standard forms."""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def __call__(self, **kwargs):
        """When this object is called, we call the default visualizing method: ``tecplot``."""
        return self.matplot(**kwargs)

    def matplot(self, density=None, i=None,
                         plot_type='contourf',
                         colormap='RdBu',
                         num_color_bar_ticks=5):
        """

        :param density:
        :param i: Plot which trace elements?
        :param plot_type: Plot type?
        :param colormap:
        :param num_color_bar_ticks:
        :return:
        """
        if density is None:
            if plot_type == 'quiver':
                density = 500
            elif plot_type == 'contourf':
                density = 10000
            else:
                raise NotImplementedError(f'3dCSCG 1Trace plot type={plot_type} is not implemented.')
        else:
            pass
        mesh = self._tf_.mesh
        density = int(np.sqrt(density/mesh.trace.elements.GLOBAL_num)) + 1
        xi = eta = sigma = np.linspace(-1, 1, density)
        xyz, v = self._tf_.reconstruct(xi, eta, sigma, i=i)
        xyz = cOmm.gather(xyz, root=mAster_rank)
        v = cOmm.gather(v, root=mAster_rank)
        if rAnk != mAster_rank: return

        XYZ = list()
        Vx = list()
        Vy = list()
        Vz = list()
        for _xyz_, _v_ in zip(xyz, v):
            for i in _xyz_:
                xyz_i = _xyz_[i]
                vx_i, vy_i, vz_i = _v_[i]

                XYZ.append(xyz_i)
                Vx.append(vx_i)
                Vy.append(vy_i)
                Vz.append(vz_i)

        Vx = np.array(Vx)
        Vy = np.array(Vy)
        Vz = np.array(Vz)
        del xyz, v

        if plot_type == 'quiver': # ================= quiver plot =====================================
            fig = plt.figure(figsize=(8, 7))

            ax = fig.add_subplot(111, projection='3d')
            for i, xyz in enumerate(XYZ):
                x, y, z = xyz
                ax.plot_surface(x, y, z, color=(0.5,0.5,0.5,0.5))
                vx = Vx[i]
                vy = Vy[i]
                vz = Vz[i]

                ax.quiver(x, y, z, vx, vy, vz, color='r', linewidth=0.5)

            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            plt.show()

        elif plot_type == 'contourf': # ================= contourf plot =====================================

            cmap = getattr(cm, colormap)

            fig = plt.figure(figsize=(15,6))

            # x-component ----------------------------------------------------------
            ax = fig.add_subplot(131, projection='3d')
            ax.view_init(45, 60)
            MAX = np.max(Vx)
            MIN = np.min(Vx)
            if MAX == MIN:
                MAX += 0.0001
            bounds = MAX - MIN
            Vx = Vx - MIN
            Vx = Vx / bounds
            ticks = np.linspace(MAX, MIN, num_color_bar_ticks)
            for i, xyz in enumerate(XYZ):
                x, y, z = xyz
                v = Vx[i]

                ax.plot_surface(x, y, z, facecolors=cmap(v))

            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array(np.array(ticks))
            cb = plt.colorbar(mappable, ax=ax, # ticks=np.linspace(0,1,num_ticks),
                              shrink=1, aspect=20,# extend='min',
                              orientation='vertical', )
            cb.ax.tick_params()  # labelsize=13.5)
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            plt.title('x-component')

            # y-component -------------------------------------------------------------
            ax = fig.add_subplot(132, projection='3d')
            ax.view_init(45, 60)
            MAX = np.max(Vy)
            MIN = np.min(Vy)
            if MAX == MIN:
                MAX += 0.0001
            bounds = MAX - MIN
            Vy = Vy - MIN
            Vy = Vy / bounds
            ticks = np.linspace(MAX, MIN, num_color_bar_ticks)
            for i, xyz in enumerate(XYZ):
                x, y, z = xyz
                v = Vy[i]

                ax.plot_surface(x, y, z, facecolors=cmap(v))

            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array(np.array(ticks))
            cb = plt.colorbar(mappable, ax=ax,  # ticks=np.linspace(0,1,num_ticks),
                              shrink=1, aspect=20,  # extend='min',
                              orientation='vertical', )
            cb.ax.tick_params()  # labelsize=13.5)
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            plt.title('y-component')

            # z-component -------------------------------------------------------
            ax = fig.add_subplot(133, projection='3d')
            ax.view_init(45, 60)
            MAX = np.max(Vz)
            MIN = np.min(Vz)
            if MAX == MIN:
                MAX += 0.0001
            bounds = MAX - MIN
            Vz = Vz - MIN
            Vz = Vz / bounds
            ticks = np.linspace(MAX, MIN, num_color_bar_ticks)
            for i, xyz in enumerate(XYZ):
                x, y, z = xyz
                v = Vz[i]

                ax.plot_surface(x, y, z, facecolors=cmap(v))

            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array(np.array(ticks))
            cb = plt.colorbar(mappable, ax=ax,  # ticks=np.linspace(0,1,num_ticks),
                              shrink=1, aspect=20,  # extend='min',
                              orientation='vertical', )
            cb.ax.tick_params()  # labelsize=13.5)
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            plt.title('z-component')

            plt.show()
        else:
            raise NotImplementedError(f'3dCSCG 1Trace plot type={plot_type} is not implemented.')