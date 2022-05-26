""""""


from screws.freeze.base import FrozenOnly



class iR_Gathering_Matrix(FrozenOnly):
    """Irregular gathering matrix: in each row, the local dofs could be different!"""
    def __init__(self):

        self._freeze_self_()


    @property
    def ___Pr_IS_regular___(self):
        return False