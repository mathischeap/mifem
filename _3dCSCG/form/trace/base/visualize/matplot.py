
import sys
if './' not in sys.path: sys.path.append('/')

from root.config import *
from SCREWS.frozen import FrozenOnly

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



class _3dCSCG_trace_form_Matplot(FrozenOnly):
    """"""

    def __init__(self, tf):
        """ """
        assert '3dCSCG_trace_form' in tf.standard_properties.tags
        self._tf_ = tf

        self._freeze_self_()

    def __call__(self, **kwargs):
        """When call it, we by default plot the whole field of the form. Tdo special plot, call the particular method."""
        getattr(self, f"_matplot_{self._tf_.k}Trace_")(**kwargs)

    def _matplot_0Trace_(self, **kwargs):
        return self._matplot_2Trace_(**kwargs)

    def _matplot_2Trace_(self, density=100000, i=None,
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

    def _matplot_1Trace_(self, density=None, i=None,
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