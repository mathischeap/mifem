# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly

import matplotlib.pyplot as plt
from matplotlib import cm
from root.config.main import *


class _3dCSCG_VectorField_matplot_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.boundary_values(*args, **kwargs)

    def boundary_values(
            self,
            density=5000, colormap='coolwarm',
            num_color_bar_ticks=5
    ):
        """"""
        mesh = self._f_.mesh
        NUM = mesh.trace.elements.global_num
        density = int((density/NUM)**0.5 + 1)

        x = y = z = np.linspace(-1, 1, density)
        xyz, v = self._f_.reconstruct(x, y, z, i='on_mesh_boundaries', where='trace-element')
        xyz = COMM.gather(xyz, root=MASTER_RANK)
        v = COMM.gather(v, root=MASTER_RANK)

        if RANK != MASTER_RANK:
            return

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
        VVV = [Vx, Vy, Vz]
        ticks = list()
        for i, V in enumerate(VVV):
            MAX = np.max(V)
            MIN = np.min(V)
            if MAX == MIN:
                MAX += 0.0001

            bounds = MAX - MIN
            V = V - MIN
            V = V / bounds

            VVV[i] = V
            ticks.append(np.linspace(MAX, MIN, num_color_bar_ticks))

        del Vx, Vy, Vz
        cmap = getattr(cm, colormap)
        fig = plt.figure(figsize=(15, 5))

        # --------- x - component --------------------------------------------------------------------------
        ax = fig.add_subplot(131, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(XYZ):
            x, y, z = xyz
            v = VVV[0][i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks[0]))
        cb = plt.colorbar(mappable, ax=ax,  # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,  # extend='min',
                          orientation='vertical', )

        cb.set_label(r'x-component')   # labelpad=10, size=15)
        cb.ax.tick_params()   # labelsize=13.5
        ax.set_xlabel(r'$x$', fontsize=10)
        ax.set_ylabel(r'$y$', fontsize=10)
        ax.set_zlabel(r'$z$', fontsize=10)

        # --------- y - component --------------------------------------------------------------------------
        ax = fig.add_subplot(132, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(XYZ):
            x, y, z = xyz
            v = VVV[1][i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks[1]))
        cb = plt.colorbar(mappable, ax=ax,  # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,  # extend='min',
                          orientation='vertical', )
        cb.set_label(r'y-component')  # labelpad=10, size=15)
        cb.ax.tick_params()   # labelsize=13.5
        ax.set_xlabel(r'$x$', fontsize=10)
        ax.set_ylabel(r'$y$', fontsize=10)
        ax.set_zlabel(r'$z$', fontsize=10)

        # --------- z - component --------------------------------------------------------------------------
        ax = fig.add_subplot(133, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(XYZ):
            x, y, z = xyz
            v = VVV[2][i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks[2]))
        cb = plt.colorbar(mappable, ax=ax,   # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,  # extend='min',
                          orientation='vertical', )
        cb.set_label(r'z-component')  # labelpad=10, size=15)
        cb.ax.tick_params()   # labelsize=13.5
        ax.set_xlabel(r'$x$', fontsize=10)
        ax.set_ylabel(r'$y$', fontsize=10)
        ax.set_zlabel(r'$z$', fontsize=10)
        # ==================================================================================================

        plt.show()
        return fig
