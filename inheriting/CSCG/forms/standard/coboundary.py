


from screws.freeze.inheriting.frozen_only import FrozenOnly


class CSCG_Standard_Form_Coboundary_BASE(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._next_form_ = None
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._incidenceMatrix_ = None

    @property
    def incidence_matrix(self):
        raise NotImplementedError()

    def ___PRIVATE_next_class___(self):
        raise NotImplementedError()

    def __call__(self):
        """
        When we call the coboundary object, we do the ``coboundary`` process; let ``self`` be a ``k``-form,
        it returns a ``(k+1)``-form.

        :return: A new standard ``(k+1)``-form.
        :raise AssertionError: If ``self`` has no cochain.
        """
        assert self._sf_.cochain.local is not None, "I need a cochain to perform coboundary."
        nextFmClass = self.___PRIVATE_next_class___()
        nextFmInstance = nextFmClass(
            self._sf_.mesh, self._sf_.space,
            is_hybrid = self._sf_.IS.hybrid,
            numbering_parameters = self._sf_.numbering._numbering_parameters_,
            name = 'd(' + self._sf_.standard_properties.name + ')'
        )
        selfCochain = self._sf_.cochain.local
        nextCochain = dict()
        incidence_matrix = self.incidence_matrix
        for i in self._sf_.mesh.elements:
            nextCochain[i] = incidence_matrix[i] @ selfCochain[i]
        nextFmInstance.cochain.local = nextCochain
        return nextFmInstance

    @property
    def cochain_local(self):
        """
        Return the local cochain (not the form) of its coboundary.

        :return:
        """
        selfCochain = self._sf_.cochain.local
        nextCochain = dict()
        incidence_matrix = self.incidence_matrix
        for i in self._sf_.mesh.elements:
            nextCochain[i] = incidence_matrix[i] @ selfCochain[i]
        return nextCochain
