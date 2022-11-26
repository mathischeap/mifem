



from components.freeze.base import FrozenOnly


from objects.CSCG.base.ADF.base.BC.interpret.main import CSCG_AFORM_BC_Interpret

class CSCG_ADForm_BC(FrozenOnly):
    """"""
    def __init__(self, adf):
        self._adf_ = adf
        self._boundaries_ = None
        self._involved_element_parts_ = None
        self._interpret_ = CSCG_AFORM_BC_Interpret(adf)
        self._freeze_self_()

    @property
    def CF(self):
        return self._adf_.prime.BC.CF


    @property
    def boundaries(self):
        """The valid boundaries of the BC of this form."""
        return self._boundaries_

    @boundaries.setter
    def boundaries(self, bns):
        """"""

        BNS = self._adf_.mesh.boundaries.names
        if isinstance(bns, str):
            bns = [bns,]
        else:
            pass
        assert isinstance(bns, (list, tuple)), f"pls put boundary names into a list or tuple."
        for i, bn in enumerate(bns):
            assert bn in BNS, f"boundary names [{i}] = {bn} is not a valid boundary name."
        self.___Pr_parse_involved_element_parts___(bns)
        self.interpret.RESET_cache() # reset current interpretations.
        self._boundaries_ = bns


    def ___Pr_parse_involved_element_parts___(self, bns):
        """"""
        mesh = self._adf_.mesh
        if mesh.ndim == 3:
            Res = mesh.boundaries.range_of_element_sides
        elif mesh.ndim == 2:
            Res = mesh.boundaries.range_of_element_edges
        else:
            raise Exception()

        self._involved_element_parts_ = list()
        for bn in bns:
            self._involved_element_parts_.extend(Res[bn])

    @property
    def interpret(self):
        return self._interpret_