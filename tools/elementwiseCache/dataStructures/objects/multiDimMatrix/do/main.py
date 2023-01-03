# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/18 22:19
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

import numpy as np
from root.config.main import COMM, MASTER_RANK, RANK

from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.do.helpers.eliminate_helper import \
    ___PRIVATE_eliminate_caller___
from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.do.helpers.sparse_matrix_helper import \
    ____PRIVATE_SparseMatrix_caller___
from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.do.helpers.evaluate_remaining_1 import nLS_DoEvaRM1
from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.do.helpers.reduce_to_vector_helper import \
    nLS_DoR2V_Helper

from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector


class MDM_Do(FrozenOnly):
    """"""
    def __init__(self, MDM):
        """"""
        self._MDM_ = MDM
        self._freeze_self_()

    def eliminate(self, dim, form):
        """We eliminate the dimension #dim of the MDM with the local cochain of form.

        Parameters
        ----------
        dim
        form

        Returns
        -------

        """

        caller = ___PRIVATE_eliminate_caller___(self._MDM_, dim, form)
        MDM_class = self._MDM_.__class__

        COR = list()
        for i, cor in enumerate(self._MDM_.correspondence):
            if i != dim:
                COR.append(cor)

        return MDM_class(self._MDM_.iterator, caller, COR, 'no_cache')

    def sparse_matrix(self, transpose=False):
        """Given that the mdm is a 2-d one, we can interpret it as a EWC sparse matrix.

        Parameters
        ----------
        transpose : bool
            Consider that the mdm is indexed as 'ij', the SparseMatrix will be 'ij' (False) or 'ji' (True)?

        Returns
        -------

        """
        assert self._MDM_.ndim == 2, f"only 2-dimensional MDM can be interpreted as a sparse matrix."

        caller = ____PRIVATE_SparseMatrix_caller___(self._MDM_, transpose)

        SpaMat = EWC_SparseMatrix(self._MDM_.iterator, caller, cache_key_generator='no_cache')

        f0, f1 = self._MDM_.correspondence

        if transpose:
            SpaMat.gathering_matrices = (f1, f0)
        else:
            SpaMat.gathering_matrices = (f0, f1)

        return SpaMat

    def derivative_contribution_to(self, test_variable, unknown_variable_pair, *known_pairs):
        """We will result in a EWC_SparseMatrix whose rows refer to test_variable, cols refer to
        unknown_variables, and all other dimensions will be eliminated according to the
        known_pairs.

        Parameters
        ----------
        test_variable
        unknown_variable_pair
        known_pairs : list
            A list of tuples, for example:
                [(u, shadow_u), (w, shadow_w), ...]
            This means u, w, ... are the forms in the correspondence, and we will use the cochains
            of shadow_u, shadow_vm ... to eliminate these dimensions.

        Returns
        -------

        """
        # ------ check the test_variable -----------------------------------------------
        assert test_variable in self._MDM_.correspondence, f"test_variable is not found."

        # --if the unknown variable is not in the correspondence, return None --
        unknown_variable = unknown_variable_pair[0]
        if unknown_variable not in self._MDM_.correspondence:
            return None
        # ------ we clear not relative pairs -----------------------------------------
        _pairs = list()
        for kp in known_pairs:
            if kp[0] in self._MDM_.correspondence:
                _pairs.append(kp)
        known_pairs = _pairs
        # --------------------------------------------------------------------------
        if self._MDM_.ndim == 3:

            if self._MDM_.whether.inhomogeneous:
                # noinspection PyTypeChecker
                assert len(known_pairs) == 1, \
                    f"for a 3d inhomogeneous MDM, please provide one known pair, " \
                    f"now it is {len(known_pairs)}."
                # noinspection PyUnresolvedReferences
                f, shadow_f = known_pairs[0]

                assert f in self._MDM_.correspondence
                index = self._MDM_.correspondence.index(f)
                MDM_2d = self.eliminate(index, shadow_f)

                if MDM_2d.correspondence == [test_variable, unknown_variable]:
                    SpaMat = MDM_2d.do.sparse_matrix(transpose=False)
                elif MDM_2d.correspondence == [unknown_variable, test_variable]:
                    SpaMat = MDM_2d.do.sparse_matrix(transpose=True)
                else:
                    raise Exception()

                return SpaMat

            else:
                raise NotImplementedError()

        else:
            raise NotImplementedError()

    def evaluate(self, known_paris):
        """

        Parameters
        ----------
        known_paris : list
            For example:
                [(u, shadow_u), (w, shadow_w), ...]
            We use u, w, ... the locate the forms, and we use cochains of
            shadow_u, shadow_w, ... to compute.

        Returns
        -------

        """
        eliminated_dimensions = dict()
        for pair in known_paris:
            u, sdu = pair
            if u in self._MDM_.correspondence:
                i = self._MDM_.correspondence.index(u)
                eliminated_dimensions[i] = sdu

        remaining_dimensions = self._MDM_.ndim - len(eliminated_dimensions)
        assert remaining_dimensions >= 0

        # ------ 0 ----------------------------------------------------------------------
        if remaining_dimensions == 0:
            # all dimensions being eliminated; return a real number
            MDM = self._MDM_
            _indices = 'abcdefghijklmnopq'[:MDM.ndim]
            # noinspection PyUnboundLocalVariable
            _eli_indices = list()
            _eli_form = list()
            for i in eliminated_dimensions:
                _eli_indices.append(_indices[i])
                _eli_form.append(eliminated_dimensions[i])
            _ein_str = _indices + ',' + ','.join(_eli_indices) + '->'

            local_value = list()
            for basic_unit in MDM:
                _ = np.einsum(
                    _ein_str,
                    MDM[basic_unit],
                    *[_.cochain.local[basic_unit] for _ in _eli_form],
                    optimize='optimal'
                )
                local_value.append(_)

            local_value = np.sum(local_value)

            local_value = COMM.gather(local_value, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                local_value = np.sum(local_value)
            else:
                pass
            local_value = COMM.bcast(local_value, root=MASTER_RANK)
            return local_value

        # ---------- 1 ---------------------------------------------------------------------------
        elif remaining_dimensions == 1:

            caller = nLS_DoEvaRM1(self._MDM_, eliminated_dimensions)
            SpaVec = EWC_ColumnVector(self._MDM_.iterator, caller, cache_key_generator='no_cache')
            return SpaVec

        # ========================================================================================
        else:
            raise NotImplementedError(f"remaining dimensions = {remaining_dimensions} not implemented.")

    def reduce_to_vector(self, v):
        """We reduce the MDM to a vector according to form `v` which must be an entry of the
        correspondence.

        We will use the cochain of other correspondences in real time.

        Parameters
        ----------
        v

        Returns
        -------

        """
        assert self._MDM_.ndim >= 2, f"MDM must be of dimensions >= 2, now it is {self._MDM_.ndim}."
        assert v in self._MDM_.correspondence, f"vector form is not found."
        real_time_caller = nLS_DoR2V_Helper(self._MDM_, v)
        return EWC_ColumnVector(self._MDM_.iterator, real_time_caller, cache_key_generator='no_cache')


if __name__ == '__main__':
    # mpiexec -n 4 python tools/elementwiseCache/dataStructures/objects/multiDimMatrix/do/main.py
    from __init__ import cscg3

    mesh = cscg3.mesh('cuboid', region_layout=[2, 2, 2])([3, 3, 3])
    space = cscg3.space('polynomials')((2, 2, 2))
    FC = cscg3.form(mesh, space)

    def u(t, x, y, z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t, x, y, z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t, x, y, z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u, v, w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    f1 = FC('1-f', hybrid=False)
    u2 = FC('2-f', hybrid=False)
    v2 = FC('2-f', hybrid=False)

    f1.CF = velocity
    f1.CF.current_time = 0
    f1.discretize()

    MDM = f1.special.cross_product_2f__ip_2f(u2, v2, output='MDM')
    mdm = MDM.do.eliminate(0, f1)

    mdm.do.sparse_matrix(transpose=True)
