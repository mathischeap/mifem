"""
A class that represent a part of the cochains of a form. Therefore, it
will have two major properties representing: (i) the dofs (numbering), (ii) the cochain of these
dofs. These two properties are named: dofs (an instance of PartialDofs) and cochain.

"""

import sys
if './' not in sys.path: sys.path.append('/')

from objects.CSCG.base.forms.base.BC.partial_cochain.partial_dofs.main import PartialDofs
from screws.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.BC.partial_cochain.include_from import _PartialCochain_Include_from_
from objects.CSCG.base.forms.base.BC.partial_cochain.interpretation import _PartialCochain_Interpretation_


class PartialCochain(FrozenOnly):
    """"""
    def __init__(self, CSCG_form):
        assert 'CSCG_form' in CSCG_form.standard_properties.tags
        self._form_ = CSCG_form
        self._mesh_ = CSCG_form.mesh
        self._dofs_ = PartialDofs(CSCG_form)
        self._cochain_ = dict()
        self._include_ = _PartialCochain_Include_from_(self)
        self._interpretation_ = _PartialCochain_Interpretation_(self)
        self._freeze_self_()

    @property
    def dofs(self):
        """The dofs are included in a PartialDofs instance."""
        return self._dofs_

    @property
    def cochain(self):
        """A dict: keys are mesh element numbers, values are the local cochains. So it is corresponding to
        `self.dofs.interpreted_as.local_dofs`.
        """
        return self._cochain_

    @property
    def include(self):
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
        """Return the "indicators" (not local dofs) and the "local cochains" of the involved mesh element #e."""
        return self._dofs_[e], self._cochain_[e]

    @property
    def interpreted_as(self):
        return self._interpretation_






if __name__ == '__main__':
    # mpiexec -n 4 python TOOLS\CSCG\partial_cochain\main.py
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



    pc = PartialCochain(t0)
    pc.include.boundaries('North')


    for e in pc:
        print(pc[e])

    # pd = PartialDofs(t2)
    # pd.include.boundaries(['Front', ])
    #
    # pd = PartialDofs(t0)
    # pd.include.boundaries(['Front', ])