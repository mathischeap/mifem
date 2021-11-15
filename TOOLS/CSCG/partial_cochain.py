"""
A class that represent a part of the cochains of a form. Therefore, it
will have two major properties representing: (i) the dofs (numbering), (ii) the cochain of these
dofs. These two properties are named: dofs (an instance of PartialDofs) and cochain.

"""

import sys
if './' not in sys.path: sys.path.append('./')

from TOOLS.CSCG.partial_dofs import PartialDofs
from SCREWS.frozen import FrozenOnly


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
        """Go through all involved local element numbers."""
        for e in self._dofs_:
            yield e

    def __len__(self):
        """How many involved local elements?"""
        return len(self._dofs_)

    def __contains__(self, e):
        """If element #e is involved?"""
        return e in self._dofs_

    def __getitem__(self, e):
        """Return the "indicators" (not local dofs) and the "local cochains" of the involved element #e."""
        return self._dofs_[e], self._cochain_[e]

    @property
    def interpreted_as(self):
        return self._interpretation_




class _PartialCochain_Include_from_(FrozenOnly):
    """"""
    def __init__(self, pc):
        self._pc_ = pc
        self._pd_ = pc.dofs
        self._cochain_ = pc._cochain_
        self._freeze_self_()

    def boundaries(self, boundary_names):
        """Include dofs (type-1 only) from boundaries named `boundary_names`."""
        pd = self._pd_
        local_dofs_indicators = pd.include.boundaries(boundary_names)
        f = self._pc_._form_

        cochain_type, cochain = f.discretize(target='BC')

        if cochain_type == 'locally full local cochain': # cochain.local, and locally full for all dofs in mesh elements.
            for i in local_dofs_indicators:
                if i not in self._cochain_:
                    self._cochain_[i] = list()
                assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"

                for side in local_dofs_indicators[i]:
                    dofs = f.numbering.DO. \
                        FIND.local_dofs_on_element_side(side)

                    self._cochain_[i].extend(cochain[i][dofs])

        elif cochain_type == 'Boundary only local cochain': # only cochains for dofs on mesh element side (boundary of the mesh).
            for i in local_dofs_indicators:
                assert i in cochain, \
                    f"element {i} in not in the local cochain, most likely," \
                    f"the boundaries in the func do not cover BC.valid_boundaries: {f.BC.valid_boundaries}."

                if i not in self._cochain_:
                    self._cochain_[i] = list()
                assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"

                for side in local_dofs_indicators[i]:
                    self._cochain_[i].extend(cochain[i][side])

        elif cochain_type == 'locally full local TEW cochain': # cochain.local_TEW, and locally full for all dofs on trace elements.
            t = f # we rename it because it must be a trace form.
            TEM = t.mesh.trace.elements.map

            if t.ndim == 2:
                side_index = {'U':0, 'D':1, 'L':2, 'R':3}
            elif t.ndim == 3:
                side_index = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}
            else:
                raise Exception()

            for i in local_dofs_indicators:
                if i not in self._cochain_:
                    self._cochain_[i] = list()

                assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"
                for side in local_dofs_indicators[i]:
                    TE = TEM[i][side_index[side]]

                    self._cochain_[i].extend(cochain[TE])

        else:
            raise NotImplementedError(f"Can not handle cochain_type={cochain_type}.")





class _PartialCochain_Interpretation_(FrozenOnly):
    """"""
    def __init__(self, pc):
        self._pc_ = pc
        self._freeze_self_()


    @property
    def local_cochains(self):
        """The cochain property actually is already interpreted as local cochains."""
        return self._pc_.cochain


if __name__ == '__main__':
    # mpiexec -n 4 python TOOLS\CSCG\partial_cochain.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector
    # from root.config import *

    mesh = MeshGenerator('crazy', c=0.0)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',3), ('Lobatto',4)])
    FC = FormCaller(mesh, space)

    # es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    f1 = FC('1-f', is_hybrid=False)
    t0 = FC('0-t')
    t2 = FC('2-t')

    # f1.TW.func.DO.set_func_body_as(es, 'velocity')
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