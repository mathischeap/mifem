"""
A class that represent a part of the dofs of a form.

"""
import sys
if './' not in sys.path: sys.path.append('/')



from screws.freeze.main import FrozenOnly
from typing import Dict
from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.include_from import _PartialDofs_Include_from_
from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.interpretation.main import _PartialDofs_Interpretation_

class PartialDofs(FrozenOnly):
    """"""
    def __init__(self, CSCG_form):
        assert 'CSCG_form' in CSCG_form.standard_properties.tags
        self._form_ = CSCG_form
        self._mesh_ = CSCG_form.mesh
        self._dofs_: Dict[int, list] = dict()
        self._include_ = _PartialDofs_Include_from_(self)
        self._interpretation_ = _PartialDofs_Interpretation_(self)
        self._freeze_self_()

    @property
    def dofs(self):
        """Return a dict. The keys are mesh element numbers. And values
        are lists of indicators indicating which local dofs.

        indicators:
            type-1: '1-' + 'N', 'S', 'W', 'E', 'B', 'F' for 3D or 'U', 'D', 'L', 'R'
                for 2D CSCG mesh --- the dofs on mesh element side.

        """
        return self._dofs_

    @property
    def include(self):
        """Methods used to include local dofs."""
        return self._include_

    def __iter__(self):
        """Go through all involved local mesh element numbers."""
        for e in self._dofs_:
            yield e

    def __len__(self):
        """How many involved local mesh elements?"""
        return len(self._dofs_)

    def __contains__(self, e):
        """If mesh element #e is involved?"""
        return e in self._dofs_

    def __getitem__(self, e):
        """Return the indicators for the involved mesh element #e."""
        return self._dofs_[e]

    @property
    def interpreted_as(self):
        return self._interpretation_





if __name__ == '__main__':
    # mpiexec -n 4 python TOOLS\CSCG\partial_cochain\partial_dofs\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector
    # from root.config import *

    mesh = MeshGenerator('crazy', c=0.0)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',3), ('Lobatto',4)])
    FC = FormCaller(mesh, space)

    # es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    f1 = FC('1-f', is_hybrid=False)
    t0 = FC('0-t')
    t2 = FC('2-t')

    # f1.TW.func.do.set_func_body_as(es, 'velocity')
    # f1.TW.current_time = 0
    # f1.TW.___DO_push_all_to_instant___()
    # f1.discretize()

    # pd = PartialDofs(f1)
    # pd.include.boundaries('North')
    #
    # pd = PartialDofs(t2)
    # pd.include.boundaries(['Front', ])

    pd = PartialDofs(t0)
    pd.include.boundaries(['Front', 'West'])

    local = pd.interpreted_as.local_dofs

    for i in local:
        print(i, local[i])