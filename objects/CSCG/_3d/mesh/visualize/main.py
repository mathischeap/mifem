
import sys
if './' not in sys.path: sys.path.append('./')



from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.visualize.matplot import _3dCSCG_Mesh_Visualize_Matplot




class _3dCSCG_Mesh_Visualize(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', " <MeshVisualize> "
        assert mesh.ndim == 3, " <MeshVisualize> "
        self._mesh_ = mesh
        self._matplot_ = _3dCSCG_Mesh_Visualize_Matplot(mesh)
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    @property
    def matplot(self):
        return self._matplot_







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\mesh\visualize\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 3, 3]
    mesh = MeshGenerator('crazy', c=0.25)(elements)
    # mesh = MeshGenerator('crazy_periodic', c=0.0, bounds=([0, 3], [0, 3], [0, 3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    # mesh = MeshGenerator('psc')(elements)
    mesh.visualize.matplot.grid()
    mesh.visualize.matplot.connection()
    mesh.visualize.matplot.element_distribution()
    mesh.domain.visualize()
    mesh.domain.regions.visualize()