
import sys
if './' not in sys.path: sys.path.append('./')


from SCREWS.frozen import FrozenOnly
from root.config import *

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class _3dCSCG_Field_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._matplot_ = _3dCSCG_Field_matplot_Visualize(f)
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_





class _3dCSCG_Field_matplot_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.boundary_values(*args, **kwargs)

    def boundary_values(self, *args, **kwargs):
        """"""
        if self._f_.__class__.__name__ == "_3dCSCG_ScalarField":
            return self.___PLOT_boundary_values_of_scalar___(*args, **kwargs)
        elif self._f_.__class__.__name__ == "_3dCSCG_VectorField":
            return self.___PLOT_boundary_values_of_vector___(*args, **kwargs)
        else:
            raise NotImplementedError(f"can not matplot boundary values of {self._f_}.")


    def ___PLOT_boundary_values_of_scalar___(self,
        density=500000, colormap='coolwarm',
        num_color_bar_ticks=5):
        """"""
        mesh = self._f_.mesh
        NUM = mesh.trace.elements.GLOBAL_num
        density = int((density/NUM)**0.5 + 1)

        x = y = z = np.linspace(-1, 1, density)

        xyz, v = self._f_.reconstruct(x, y, z, where='trace-element')

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
        cb.ax.tick_params(labelsize=13.5)


        plt.show()



    def ___PLOT_boundary_values_of_vector___(self,
        density=500000, colormap='coolwarm',
        num_color_bar_ticks=5):
        """"""
        mesh = self._f_.mesh
        NUM = mesh.trace.elements.GLOBAL_num
        density = int((density/NUM)**0.5 + 1)

        x = y = z = np.linspace(-1, 1, density)
        xyz, v = self._f_.reconstruct(x, y, z, where='trace-element')
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
        fig = plt.figure(figsize=(15,5))

        # --------- x - component --------------------------------------------------------------------------
        ax = fig.add_subplot(131, projection='3d')
        ax.view_init(45, 60)

        for i, xyz in enumerate(XYZ):
            x, y, z = xyz
            v = VVV[0][i]

            ax.plot_surface(x, y, z, facecolors=cmap(v))

        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array(np.array(ticks[0]))
        cb = plt.colorbar(mappable, ax=ax, # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,# extend='min',
                          orientation='vertical', )
        cb.set_label(r'x-component')#,labelpad=10, size=15)
        cb.ax.tick_params() # labelsize=13.5
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
        cb = plt.colorbar(mappable, ax=ax, # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,# extend='min',
                          orientation='vertical', )
        cb.set_label(r'y-component')#,labelpad=10, size=15)
        cb.ax.tick_params() #  labelsize=13.5
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
        cb = plt.colorbar(mappable, ax=ax, # ticks=np.linspace(0,1,num_ticks),
                          shrink=1, aspect=20,# extend='min',
                          orientation='vertical', )
        cb.set_label(r'z-component')#,labelpad=10, size=15)
        cb.ax.tick_params() # labelsize=13.5
        ax.set_xlabel(r'$x$', fontsize=10)
        ax.set_ylabel(r'$y$', fontsize=10)
        ax.set_zlabel(r'$z$', fontsize=10)
        #===================================================================================================

        plt.show()
        return fig


if __name__ == '__main__':

    # mpiexec -n 6 python _3dCSCG\field\visualize.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)
    SS = FC('scalar', 1)
    BS = FC('scalar', {'North': 1, 'West':p})

    BS.current_time=0
    BS.visualize()
