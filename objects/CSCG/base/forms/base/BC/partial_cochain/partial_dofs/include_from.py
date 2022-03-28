




from screws.freeze.main import FrozenOnly


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