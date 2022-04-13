


from screws.freeze.base import FrozenOnly




class SpaMat_Adjust(FrozenOnly):
    """Like Customize, Adjust will make the changes in real time and make a new EWC_SparseMatrix
    immediately.
    """
    def __init__(self, mat):
        """"""
        self._mat_ = mat
        self._freeze_self_()

