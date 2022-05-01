


from screws.freeze.base import FrozenOnly
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate



class EWC_ColVec_Blocks(FrozenOnly):
    """"""
    def __init__(self, col_vec):
        if col_vec.con_shape is False:
            raise Exception(f"Blocks only valid for concatenate EWC_ColumnVector.")
        assert col_vec.gathering_matrix is not None, f"To access blocks, first set GM."
        self._shape_ = col_vec.con_shape
        self.___BLOCKS___ = col_vec._DG_.vectors
        self._col_vec_ = col_vec
        self._freeze_self_()


    @property
    def shape(self):
        return self._shape_

    def __getitem__(self, item):
        B = self.___BLOCKS___[item]
        return concatenate(B)
