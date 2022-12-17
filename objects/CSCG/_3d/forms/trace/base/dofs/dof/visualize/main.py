
import sys
if './' not in sys.path:
    sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.trace.base.dofs.dof.visualize.matplot._0tr import _3dCSCG_0TF_DOF_Matplot
from objects.CSCG._3d.forms.trace.base.dofs.dof.visualize.matplot._1tr import _3dCSCG_1TF_DOF_Matplot
from objects.CSCG._3d.forms.trace.base.dofs.dof.visualize.matplot._2tr import _3dCSCG_2TF_DOF_Matplot


class _3dCSCG_Trace_forms_DOF_VISUALIZE(FrozenOnly):
    """A dof of a trace form."""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._tf_ = dof._tf_
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """Call the default plotter: The matplot."""
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        """The plotter uses the matplotlib library."""
        if self._matplot_ is None:

            k = self._tf_.k
            if k == 0:
                self._matplot_ = _3dCSCG_0TF_DOF_Matplot(self._dof_)
            elif k == 0:
                self._matplot_ = _3dCSCG_1TF_DOF_Matplot(self._dof_)
            elif k == 0:
                self._matplot_ = _3dCSCG_2TF_DOF_Matplot(self._dof_)
            else:
                raise Exception()

        return self._matplot_


if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\trace\dofs\dof\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.3)(None)
    space = SpaceInvoker('polynomials')([('Lobatto', 1), ('Lobatto', 1), ('Lobatto', 1)])
    FC = FormCaller(mesh, space)


    t0 = FC('0-t')

    dofs = t0.dofs
    DI = dofs[0]
    bf = DI.basis_function

    vis = DI.visualize()
