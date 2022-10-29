# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:14 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard._2.base.main import miUsTriangular_S2F_Base


class miUsTriangular_S2F_Inner(miUsTriangular_S2F_Base):
    """"""

    def __init__(self, mesh, space, name='Tri-is2f'):
        """"""
        super(miUsTriangular_S2F_Inner, self).__init__(mesh, space, 'inner', name)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_2/inner/main.py
    import numpy as np
    from root.config.main import RANK
    import matplotlib.pyplot as plt
    from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
    from objects.miUsGrid.triangular.__test__.Random.test_mesh import mesh
    from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

    P = [i for i in range(1,30)]
    Error = list()
    for p in P:

        space = miUsGrid_TriangularFunctionSpace(p)

        f2 = miUsTriangular_S2F_Inner(mesh, space)

        def func(t, x, y): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + t

        scalar = miUsGrid_Triangular_Scalar(mesh, func)

        f2.CF = scalar

        scalar.current_time = 0

        f2.discretize()

        error = f2.error.L()

        Error.append(error)

    if RANK == 0:
        plt.rc('text', usetex=True)
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
        figure = plt.figure(figsize=(5,4))
        plt.semilogy(P, Error, '-k', linewidth=1)
        plt.xlabel('N')
        plt.ylabel(r"$\left\|f-f^h\right\|_{L^2}$")
        plt.show()