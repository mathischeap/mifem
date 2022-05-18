# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 7:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly




class _2nCSCG_RF2_VectorVisualize(FrozenOnly):
    """"""

    def __init__(self, cf):
        """"""
        self._cf_ = cf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)


    def matplot(self, *args, **kwargs):
        """"""
        ftype = self._cf_.ftype

        if ftype == 'standard':
            return self.___Pr_matplot_standard___(*args, **kwargs)
        else:
            raise NotImplementedError(f"{ftype}.")


    def ___Pr_matplot_standard___(self, density=None, plot_type='contourf', **kwargs):
        """

        Parameters
        ----------
        density
        kwargs

        Returns
        -------

        """
        if density is None:
            if plot_type in ('contourf', 'contour'):
                density = 20
            elif plot_type == 'quiver':
                density = 5
            else:
                density=15
        else:
            pass
        mesh = self._cf_.mesh
        coo = mesh.coordinates.homogeneous(density, ndim=2)
        xy, v = self._cf_.reconstruct(coo, ravel=False)
        v.visualize(xy, plot_type=plot_type, **kwargs)







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/fields/vector/visualize.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    from objects.nCSCG.rf2._2d.fields.vector.main import _2nCSCG_RF2_VectorField

    mesh = rm2(100)

    import numpy as np
    def p(t, x, y): return 2 * np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return 2 * np.cos(np.pi*x) * np.sin(np.pi*y) + t

    v = _2nCSCG_RF2_VectorField(mesh, (p,q))
    v.current_time = 0

    v.visualize(plot_type='quiver')