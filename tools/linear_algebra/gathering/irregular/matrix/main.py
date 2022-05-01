""""""


from screws.freeze.base import FrozenOnly



class IrregularGatheringMatrix(FrozenOnly):
    """Irregular gathering matrix: in each row, the local dofs could be different!"""
    def __init__(self):

        self._freeze_self_()