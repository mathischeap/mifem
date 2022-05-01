


import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard.base.dofs.dof.visualize.matplot._0sf import \
    _3dCSCG_SF_DOF_VISUALIZE_matplot_0SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.visualize.matplot._1sf import \
    _3dCSCG_SF_DOF_VISUALIZE_matplot_1SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.visualize.matplot._2sf import \
    _3dCSCG_SF_DOF_VISUALIZE_matplot_2SF
from objects.CSCG._3d.forms.standard.base.dofs.dof.visualize.matplot._3sf import \
    _3dCSCG_SF_DOF_VISUALIZE_matplot_3SF




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

            if self._dof_._sf_.k == 0:
                self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot_0SF(self._dof_)
            elif self._dof_._sf_.k == 1:
                self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot_1SF(self._dof_)
            elif self._dof_._sf_.k == 2:
                self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot_2SF(self._dof_)
            elif self._dof_._sf_.k == 3:
                self._matplot_ = _3dCSCG_SF_DOF_VISUALIZE_matplot_3SF(self._dof_)
            else:
                raise Exception()

        return self._matplot_







if __name__ == '__main__':
    # mpiexec -n 6 python objects\CSCG\_3d\forms\standard\base\dofs\dof\visualize\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',2), ('Lobatto',1)])
    FC = FormCaller(mesh, space)

    f = FC('0-f', is_hybrid=True)

    dofs = f.dofs
    for i in dofs:
        DI = dofs[i]
        DI.visualize()