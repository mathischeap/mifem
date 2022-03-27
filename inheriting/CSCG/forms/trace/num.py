
from screws.freeze.base import FrozenOnly


class CSCG_trace_form_NUM(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    @property
    def basis(self):
        """Return number of dofs in a mesh element."""
        return self._tf_._NUM_basis_

    @property
    def basis_components(self):
        """
        (Dict[str, Tuple[int]]) Return a dict whose keys are 'N', 'S', 'W', 'E', 'B','F' and values are
        tuples of basis functions on the side. If we have a 0- or 2-trace-form, it is like
        t.NUM_basis_components['N'] = (x,) and t.NUM_basis_onside['N'] = x. While when we have a
        1-trace-form, it is like t.NUM_basis_components['N'] = (x, y) and t.NUM_basis_onside['N'] = x + y.
        """
        return self._tf_._NUM_basis_components_

    @property
    def basis_onside(self):
        """(Dict[str, int]) Return a dict. See docstring of ``NUM_basis_components`` method."""
        return self._tf_._NUM_basis_onside_