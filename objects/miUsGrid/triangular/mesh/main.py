# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.base.mesh.main import miUsGrid_MeshBase
from objects.miUsGrid.triangular.mesh.construct.main import miUsGrid_TriangularMesh_Construct
from objects.miUsGrid.triangular.mesh.elements.main import miUsGrid_TriangularMesh_Elements
from objects.miUsGrid.triangular.mesh.visualize.main import miUsGrid_TriangularMesh_Visualize
from objects.miUsGrid.triangular.mesh.boundaries.main import miUsGrid_TriangularMesh_Boundaries

class miUsGrid_TriangularMesh(miUsGrid_MeshBase):
    """"""

    def __init__(self, source, boundaries, name='NoNameMesh'):
        """

        Parameters
        ----------
        source :
            We construct a miUsGrid_TriangularMesh from this source.
        boundaries : dict
            A dict whose keys are boundary names and values are functions take coordinates
            as inputs and return bool which represent if the coordinates are on the boundaries.

            For example,

                boundaries = {
                    'Upper': fu(x, y): return x == 0,
                    'Down': fd(x, y): return x == 1,
                    'Left': fl(x, y): return y == 0,
                    'right': fr(x, y): return y == 1

                }

            This `boundaries` defines the four boundaries of the domain [0,1]^2.

        """
        super(miUsGrid_TriangularMesh, self).__init__(2, name)
        self._elements_ = miUsGrid_TriangularMesh_Elements(self)
        miUsGrid_TriangularMesh_Construct(source)(self)
        self._visualize_ = miUsGrid_TriangularMesh_Visualize(self)
        self._boundaries_ = miUsGrid_TriangularMesh_Boundaries(self, boundaries)





if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/main.py
    from objects.miUsGrid.triangular.__test__.Random.test_mesh import mesh
    # mesh.visualize()
    print(mesh == mesh)