

import numpy as np
from screws.freeze.main import FrozenOnly
from screws.exceptions import ElementEdgePairError
from objects.CSCG._2d.mesh.do.find import _2dCSCG_Mesh_DO_FIND
import random


class _2dCSCG_Mesh_DO(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._FIND_ = _2dCSCG_Mesh_DO_FIND(self)
        self._freeze_self_()

    def reset_cache(self):
        self._mesh_.___PRIVATE_reset_cache___()

    @staticmethod
    def parse_element_edge_pair(eP: str):
        """Element side pairs are also used for trace element keys."""
        if eP.count('-') == 1:
            # must be regular pair to domain boundary!
            elementOne = int(eP.split('-')[0])
            edgeOne = eP.split('-')[1][0]
            boundaryName = eP.split('|')[1]
            return 'regular|domainBoundary', elementOne, edgeOne, boundaryName
        elif eP.count('-') == 2:
            elementOne, pairTypeINFO, elementTwo = eP.split('-')
            elementOne = int(elementOne)
            elementTwo = int(elementTwo)
            if len(pairTypeINFO) == 3 and pairTypeINFO[1] == '|':
                # regular pair; conforming pair; N|S, W|E, B|F and no twist!
                edgeOne = pairTypeINFO[0]
                edgeTwo = pairTypeINFO[2]
                return 'regular|regular', \
                       [elementOne, elementTwo], edgeOne + edgeTwo, None  # None for future extension.
                # for all kinds of pair, return has followed this same rule!
            else:
                raise ElementEdgePairError(f"Pair: {pairTypeINFO} is not understandable.")
        else:
            raise Exception('elementSidePair format wrong!')



    @property
    def find(self):
        return self._FIND_



    def regionwsie_stack(self, *ndas):
        """
        We use this method to stack a ndarray regions-wise. This function is very useful
        in plotting reconstruction data. Since in a region, the date are structured mesh element wise,
        so we can only plot element by element. But if we group data from elements of the same
        region, then we can plot region by region. This very much increases the plotting speed.

        We better first collect all data to one core. That means we should always call this method
        from a single core, mostly, the master core.

        Parameters
        ----------
        ndas : ndarray
            The ndarray to be stacked. The ndim of the 'nda' must be self.ndim + 1. and
            `np.shape('nda')[0]` must == `self.elements.num`.

        Returns
        -------
        output : tuple

        """
        _SD_ = tuple()
        for nda in ndas:
            if isinstance(nda, dict):
                for _ in nda: assert np.ndim(nda[_]) == 2
            else:
                assert np.ndim(nda) == 2 + 1
            if isinstance(nda, dict):
                assert len(nda) == self._mesh_._num_total_elements_
            else:
                assert len(nda) == self._mesh_._num_total_elements_
            _sd_ = {}
            if isinstance(nda, dict):
                ij = np.shape(nda[0])
            else:
                ij = np.shape(nda)[1:]
            I, J = ij
            ALL_element_global_numbering_ = \
                self._mesh_.___PRIVATE_generate_ALL_element_global_numbering___()
            for Rn in ALL_element_global_numbering_:
                region_data_shape = [ij[i] * self._mesh_._element_layout_[Rn][i] for i in range(2)]
                _sd_[Rn] = np.zeros(region_data_shape)
                for j in range(self._mesh_._element_layout_[Rn][1]):
                    for i in range(self._mesh_._element_layout_[Rn][0]):
                        _sd_[Rn][i * I:(i + 1) * I, j * J:(j + 1) * J] = \
                            nda[ALL_element_global_numbering_[Rn][i, j]]
            _SD_ += (_sd_,)
        _SD_ = _SD_[0] if len(ndas) == 1 else _SD_
        return _SD_


    def generate_random_coordinates(self):
        """ We will generate some random coordinates with in this domain in a format of local mesh
        element wise. They are mainly for testing purposes.

        :return: two 1d array representing x, y coordinates.
        """
        amount = random.randint(10, 100) # 3*amount points in #amount local mesh elements.

        assert isinstance(amount, int) and amount > 0, \
            f"amount={amount} ({amount.__class__.__name__}) is wrong."

        mesh = self._mesh_
        num_local_mesh_elements = mesh.elements.num

        if num_local_mesh_elements == 0:
            return np.array([]), np.array([])

        elif num_local_mesh_elements <= amount: # we will use all local mesh elements
            USED = mesh.elements.indices

        else:
            USED = random.sample(mesh.elements.indices, amount)

        xi, et = [np.random.rand(3) * 2 - 1 for _ in range(2)]
        X = list()
        Y = list()
        for i in USED:
            element = mesh.elements[i]
            x, y = element.coordinate_transformation.mapping(xi, et)
            X.append(x)
            Y.append(y)
        X = np.concatenate(X)
        Y = np.concatenate(Y)

        return X, Y