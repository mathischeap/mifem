
import sys
if './' not in sys.path: sys.path.append('/')


from screws.freeze.main import FrozenOnly
from root.config.main import *
from objects.CSCG._3d.fields.vector.visualize.matplot import _3dCSCG_VectorField_matplot_Visualize






class _3dCSCG_VectorField_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._matplot_ = _3dCSCG_VectorField_matplot_Visualize(f)
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_












if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\fields\base\visualize.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)
    SS = FC('scalar', 1)
    BS = FC('scalar', {'North': 1, 'West':p})

    BS.current_time=0
    BS.visualize()
