


import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from _3dCSCG.forms.standard.base.dofs.dof.visualize.matplot import _3dCSCG_SF_DOF_VISUALIZE_matplot




class _3dCSCG_SF_DOF_VISUALIZE(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._mesh_ = dof._sf_.mesh
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot(self._dof_)
        return self._matplot_







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\standard\base\dofs\dof\visualize\main.py
    from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    f = FC('3-f', is_hybrid=True)

    # SPL = f0.numbering.sharing_physical_locations
    #
    # print(SPL)

    dofs = f.dofs
    for i in dofs:
        DI = dofs[i]
        DI.visualize()