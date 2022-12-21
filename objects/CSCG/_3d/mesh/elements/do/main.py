# -*- coding: utf-8 -*-
""""""
from components.freeze.main import FrozenOnly
from root.config.main import COMM, np
import matplotlib.pyplot as plt

from objects.CSCG._3d.mesh.elements.do.find import _3dCSCG_Mesh_Elements_DO_FIND


class _3dCSCG_Mesh_Elements_DO(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._FIND_ = None
        self._freeze_self_()

    def illustrate_element(self, i, density_factor=2):
        """We use this method to illustrate a mesh element.

        :param int i: We illustrate mesh element #i.
        :param int density_factor: How refined the plots are? be in {1,2,3,...}
        :return:

        """
        if i not in self._elements_:
            COMM.barrier()
            return

        density = 5 + 4 * density_factor
        i0 = 1 + density_factor
        i1 = 2 * density_factor + 2
        i2 = 3 * density_factor + 3
        _ = np.linspace(-1, 1, density)
        r, s = np.meshgrid(_, _, indexing='ij')
        anchors = (  # the points we will plot the outward unit norm vector
            [i0, i0],
            [i0, i2],
            [i2, i0],
            [i2, i2],
            [i1, i1],
        )
        uv_r = np.array([r[indices[0], indices[1]] for indices in anchors])
        uv_s = np.array([s[indices[0], indices[1]] for indices in anchors])

        element = self._elements_[i]
        sides = element.sides

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        x_lim, y_lim, z_lim = [list(), list(), list()]
        for sn in sides:
            side = sides[sn]
            CT = side.coordinate_transformation
            x, y, z = CT.mapping(r, s)
            ax.plot_surface(x, y, z)
            x_lim.append(np.min(x))
            x_lim.append(np.max(x))
            y_lim.append(np.min(y))
            y_lim.append(np.max(y))
            z_lim.append(np.min(z))
            z_lim.append(np.max(z))

        x_range = np.max(x_lim) - np.min(x_lim)
        y_range = np.max(y_lim) - np.min(y_lim)
        z_range = np.max(z_lim) - np.min(z_lim)

        mean_range = (x_range + y_range + z_range) / 6

        for sn in sides:
            side = sides[sn]
            CT = side.coordinate_transformation
            x, y, z = CT.mapping(uv_r, uv_s)
            u, v, w = CT.outward_unit_normal_vector(uv_r, uv_s)
            ax.quiver(x, y, z, u*mean_range, v*mean_range, w*mean_range, color='red', linewidth=0.8)

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')
        plt.title(f"mesh element: {i}")

        plt.show()
        plt.close()
        COMM.barrier()

    @property
    def find(self):
        if self._FIND_ is None:
            self._FIND_ = _3dCSCG_Mesh_Elements_DO_FIND(self._elements_)
        return self._FIND_
