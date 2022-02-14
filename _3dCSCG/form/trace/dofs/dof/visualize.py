
import sys
if './' not in sys.path: sys.path.append('./')

from SCREWS.frozen import FrozenOnly



class _3dCSCG_Trace_forms_DOF_VISUALIZE(FrozenOnly):
    """A dof of a trace form."""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._tf_ = dof._tf_
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """Call the default ploter: The matplot."""
        return self.matplot(*args, **kwargs)


    def matplot(self,  *args, **kwargs):
        """The ploter uses the matplotlib library."""
        return getattr(self, f"___PRIVATE_matplot_{self._tf_.k}trace___")(*args, **kwargs)


    def ___PRIVATE_matplot_0trace___(self):
        """"""
        print("___PRIVATE_matplot_0trace___")





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\trace\dofs\dof\visualize.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)(None)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)


    t0 = FC('0-t')

    dofs = t0.dofs
    DI = dofs[0]
    bf = DI.basis_function

    vis = DI.visualize()