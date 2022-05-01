
from screws.freeze.main import FrozenOnly




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


        if cochain_type == 'locally full local cochain':
            # cochain.local, and locally full for all dofs in mesh elements.
            if f.ndim == 3:
                for i in local_dofs_indicators:
                    if i not in self._cochain_:
                        self._cochain_[i] = list()
                    assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"

                    for side in local_dofs_indicators[i]:
                        dofs = f.numbering.do. \
                            find.local_dofs_on_element_side(side)

                        self._cochain_[i].extend(cochain[i][dofs])
            elif f.ndim == 2:
                for i in local_dofs_indicators:
                    if i not in self._cochain_:
                        self._cochain_[i] = list()
                    assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"

                    for side in local_dofs_indicators[i]:
                        dofs = f.numbering.do. \
                            find.local_dofs_on_element_edge(side)

                        self._cochain_[i].extend(cochain[i][dofs])
            else:
                raise Exception()

        elif cochain_type == 'Boundary only local cochain':
            # only cochains for dofs on mesh element side (boundary of the mesh).
            for i in local_dofs_indicators:
                assert i in cochain, \
                    f"element {i} in not in the local cochain, most likely," \
                    f"the boundaries in the func do not cover BC.valid_boundaries: {f.BC.valid_boundaries}."

                if i not in self._cochain_:
                    self._cochain_[i] = list()

                assert len(local_dofs_indicators[i]) > 0, f"empty for element #{i}"

                for side in local_dofs_indicators[i]:
                    self._cochain_[i].extend(cochain[i][side])

        elif cochain_type == 'locally full local TEW cochain':
            # cochain.local_TEW, and locally full for all dofs on trace elements.
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

