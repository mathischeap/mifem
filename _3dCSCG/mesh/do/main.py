



from screws.freeze.main import FrozenOnly
from _3dCSCG.mesh.do.find import _3dCSCG_Mesh_DO_FIND
import numpy as np





class _3dCSCG_Mesh_DO(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._FIND_ = _3dCSCG_Mesh_DO_FIND(self)
        self._freeze_self_()

    def reset_cache(self):
        self._mesh_.___PRIVATE_reset_cache___()

    def parse_element_side_pair(self, eP):
        return self._mesh_.___PRIVATE_do_parse_element_side_pair___(eP)

    @property
    def find(self):
        return self._FIND_


    def regionwsie_stack(self, *nda_s):
        """
        Wo should only use it in one core (first collect all data to this core).

        We use this method to stack a ndarray regions-wise. This function is very useful
        in plotting reconstruction data. Since in a regions, the elements are structure,
        we can plot element by element. But if we group data from elements of the same
        regions, then we can plot regions by regions. This very increase the plotting speed
        significantly.

        Parameters
        ----------
        nda_s : ndarray
            The ndarray to be stacked. The ndim of the 'nda' must be self.ndim + 1. and
            np.shape('nda')[0] must == self._num_total_elements_.

        Returns
        -------
        output : tuple
        """
        _SD_ = tuple()
        for nda in nda_s:
            assert np.ndim(nda) == self._mesh_.ndim + 1
            assert np.shape(nda)[0] == self._mesh_._num_total_elements_
            _sd_ = dict()
            ijk = np.shape(nda)[1:]
            I, J, K = ijk
            ALL_element_global_numbering_ = \
                self._mesh_.___PRIVATE_generate_ALL_element_global_numbering___()
            for Rn in ALL_element_global_numbering_:
                region_data_shape = [ijk[i] * self._mesh_._element_layout_[Rn][i] for i in range(3)]
                _sd_[Rn] = np.zeros(region_data_shape)
                for k in range(self._mesh_._element_layout_[Rn][2]):
                    for j in range(self._mesh_._element_layout_[Rn][1]):
                        for i in range(self._mesh_._element_layout_[Rn][0]):
                            _sd_[Rn][i * I:(i + 1) * I, j * J:(j + 1) * J, k * K:(k + 1) * K] = \
                                nda[ALL_element_global_numbering_[Rn][i, j, k]]
            _SD_ += (_sd_,)
        _SD_ = _SD_[0] if len(nda_s) == 1 else _SD_
        return _SD_