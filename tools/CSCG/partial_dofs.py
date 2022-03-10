"""
A class that represent a part of the dofs of a form.

"""
import sys
if './' not in sys.path: sys.path.append('./')



from screws.freeze.main import FrozenOnly
from typing import Dict


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



class _PartialDofs_Include_from_(FrozenOnly):
    """A wrapper of local dof including methods."""
    def __init__(self, pd):
        self._pd_ = pd
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._freeze_self_()

    def boundaries(self, boundary_names):
        """Include dofs (type-1 only) from boundaries named `boundary_names`."""

        if isinstance(boundary_names, str):
            boundary_names = [boundary_names,]

        mesh = self._mesh_

        for bn in boundary_names:
            assert bn in mesh.boundaries.names, \
                f"boundary named {bn} is not found!"

        Res = mesh.boundaries.range_of_element_sides

        new_added = dict()

        for bn in boundary_names:
            local_elements_and_sides = Res[bn]
            for local_element_and_side in local_elements_and_sides:
                element, side = int(local_element_and_side[:-1]), \
                                local_element_and_side[-1]
                if element not in self._pd_._dofs_:
                    self._pd_._dofs_[element] = list()
                self._pd_._dofs_[element].append('1-'+side) # type-1 indicators

                if element not in new_added:
                    new_added[element] = list()
                new_added[element].append(side) # must be type-1, we ignore '1-'.

        return new_added







class _PartialDofs_Interpretation_(FrozenOnly):
    """A wrapper of local dof interpretation methods."""
    def __init__(self, pd):
        self._pd_ = pd
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._local_dofs_ = None
        self._globe_dofs_ = None
        self._freeze_self_()

    @property
    def local_dofs(self):
        if self._local_dofs_ is None:
            self._local_dofs_ = _PartialDofs_Interpretation_Local_(self)
        return self._local_dofs_

    @property
    def globe_dofs(self):
        if self._globe_dofs_ is None:
            self._globe_dofs_ = _PartialDofs_Interpretation_Globe_(self)
        return self._globe_dofs_



class _PartialDofs_Interpretation_Local_(FrozenOnly):
    """The class of the local dof interpretation."""
    def __init__(self, interpretation):
        self._interpretation_ = interpretation
        pd = interpretation._pd_
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._dofs_ = pd._dofs_
        self._freeze_self_()

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
        """interpret the indicators for the involved element #e."""
        indicators = self._dofs_[e]
        DOFs = list()
        for indi in indicators:
            if indi[:2] == '1-': # type-1 indicators
                side = indi[2]
                # so the numbering property must have this method implemented.
                dofs = self._form_.numbering.do.\
                    find.local_dofs_on_element_side(side)
                DOFs.extend(dofs)
            else:
                raise NotImplementedError(
                    f"Cannot interpret as local_dofs for indicator "
                    f"{indi} for element #{e}.")
        return DOFs



class _PartialDofs_Interpretation_Globe_(FrozenOnly):
    """The class of the local dof interpretation."""
    def __init__(self, interpretation):
        self._interpretation_ = interpretation
        pd = interpretation._pd_
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._dofs_ = pd._dofs_
        self._freeze_self_()


    def __iter__(self):
        """Go through all involved local element numbers."""
        for e in self._dofs_:
            yield e

    def __getitem__(self, e):
        """interpret the indicators for the involved element #e."""
        raise NotImplementedError()








if __name__ == '__main__':
    # mpiexec -n 4 python TOOLS\CSCG\partial_dofs.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector
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