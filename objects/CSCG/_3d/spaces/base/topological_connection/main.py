# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/4/2022 8:28 PM
"""
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from scipy.sparse import csr_matrix
import numpy as np

class _3dCSCG_SpaceTopologicalConnections(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._mesh_ = None
        self._freeze_self_()

    def __call__(self, *forms):
        """Here we return the topological connection matrix between forms.

        Parameters
        ----------
        forms

        Returns
        -------

        """
        if len(forms) == 3:
            return self.___Pr3TC_Three___(*forms)
        else:
            raise NotImplementedError()

    def ___Pr3TC_Three___(self, *forms):
        """"""
        cns = [_.__class__.__name__ for _ in forms]
        if set(cns) == {'_3dCSCG_0LocalTrace', '_3dCSCG_0Trace', '_3dCSCG_0Form'}:
            inputs = [None, None, None]
            for i, cn in enumerate(cns):
                if cn == '_3dCSCG_0Form':
                    inputs[0] = forms[i]
                elif cn == '_3dCSCG_0LocalTrace':
                    inputs[1] = forms[i]
                elif cn == '_3dCSCG_0Trace':
                    inputs[2] = forms[i]
                else:
                    raise Exception()
            return self.___Pr3TC_0sf_0ltf_0tf___(*inputs)

        elif set(cns) == {'_3dCSCG_2LocalTrace', '_3dCSCG_2Trace', '_3dCSCG_2Form'}:
            inputs = [None, None, None]
            for i, cn in enumerate(cns):
                if cn == '_3dCSCG_2Form':
                    inputs[0] = forms[i]
                elif cn == '_3dCSCG_2LocalTrace':
                    inputs[1] = forms[i]
                elif cn == '_3dCSCG_2Trace':
                    inputs[2] = forms[i]
                else:
                    raise Exception()
            return self.___Pr3TC_2sf_2ltf_2tf___(*inputs)

        else:
            raise NotImplementedError()

    @staticmethod
    def ___Pr3TC_0sf_0ltf_0tf___(sf0, ltf0, tf0):
        """"""
        assert sf0.space == ltf0.space == tf0.space and \
               sf0.mesh == ltf0.mesh == tf0.mesh, f"mesh or space does not match."
        assert sf0.whether.hybrid and not ltf0.whether.hybrid and not tf0.whether.hybrid, \
            f"hybrid setting wrong."

        nb_sf0 = sf0.num.basis
        bn_ltf0 = ltf0.num.basis
        nb_tf0 = tf0.num.basis

        T_ltf_sf = np.zeros((bn_ltf0, nb_sf0), dtype=int)
        T_ltf_tf = np.zeros((bn_ltf0, nb_tf0), dtype=int)

        sf_local_numbering = sf0.numbering.local[0]
        ltf_lg = ltf0.numbering.local_gathering
        tf_lg = tf0.numbering.local_gathering

        T_ltf_sf[ltf_lg['N'], sf_local_numbering[0, :, :].ravel('F')] = -1
        T_ltf_sf[ltf_lg['S'], sf_local_numbering[-1, :, :].ravel('F')] = 1
        T_ltf_sf[ltf_lg['W'], sf_local_numbering[:, 0, :].ravel('F')] = -1
        T_ltf_sf[ltf_lg['E'], sf_local_numbering[:, -1, :].ravel('F')] = 1
        T_ltf_sf[ltf_lg['B'], sf_local_numbering[:, :, 0].ravel('F')] = -1
        T_ltf_sf[ltf_lg['F'], sf_local_numbering[:, :, -1].ravel('F')] = 1

        T_ltf_tf[ltf_lg['N'], tf_lg['N']] = -1
        T_ltf_tf[ltf_lg['S'], tf_lg['S']] = 1
        T_ltf_tf[ltf_lg['W'], tf_lg['W']] = -1
        T_ltf_tf[ltf_lg['E'], tf_lg['E']] = 1
        T_ltf_tf[ltf_lg['B'], tf_lg['B']] = -1
        T_ltf_tf[ltf_lg['F'], tf_lg['F']] = 1

        # same + and -, Do T_ltf_sf & -T_ltf_tf apply the correct coupling.

        return (
            EWC_SparseMatrix(sf0.mesh, csr_matrix(T_ltf_sf), 'constant'),
            EWC_SparseMatrix(sf0.mesh, csr_matrix(T_ltf_tf), 'constant')
        )

    @staticmethod
    def ___Pr3TC_2sf_2ltf_2tf___(sf2, ltf2, tf2):
        """"""
        assert sf2.space == ltf2.space == tf2.space and \
               sf2.mesh == ltf2.mesh == tf2.mesh, f"mesh or space does not match."
        assert sf2.whether.hybrid and not ltf2.whether.hybrid and not tf2.whether.hybrid, \
            f"hybrid setting wrong."

        nb_sf2 = sf2.num.basis
        bn_ltf2 = ltf2.num.basis
        nb_tf2 = tf2.num.basis

        T_ltf_sf = np.zeros((bn_ltf2, nb_sf2), dtype=int)
        T_ltf_tf = np.zeros((bn_ltf2, nb_tf2), dtype=int)

        sf_local_numbering = sf2.numbering.local
        ltf_lg = ltf2.numbering.local_gathering
        tf_lg = tf2.numbering.local_gathering

        T_ltf_sf[ltf_lg['N'], sf_local_numbering[0][0, :, :].ravel('F')] = -1
        T_ltf_sf[ltf_lg['S'], sf_local_numbering[0][-1, :, :].ravel('F')] = 1
        T_ltf_sf[ltf_lg['W'], sf_local_numbering[1][:, 0, :].ravel('F')] = -1
        T_ltf_sf[ltf_lg['E'], sf_local_numbering[1][:, -1, :].ravel('F')] = 1
        T_ltf_sf[ltf_lg['B'], sf_local_numbering[2][:, :, 0].ravel('F')] = -1
        T_ltf_sf[ltf_lg['F'], sf_local_numbering[2][:, :, -1].ravel('F')] = 1

        T_ltf_tf[ltf_lg['N'], tf_lg['N']] = -1
        T_ltf_tf[ltf_lg['S'], tf_lg['S']] = 1
        T_ltf_tf[ltf_lg['W'], tf_lg['W']] = -1
        T_ltf_tf[ltf_lg['E'], tf_lg['E']] = 1
        T_ltf_tf[ltf_lg['B'], tf_lg['B']] = -1
        T_ltf_tf[ltf_lg['F'], tf_lg['F']] = 1

        # same + and -, Do T_ltf_sf & -T_ltf_tf apply the correct coupling.

        return (
            EWC_SparseMatrix(sf2.mesh, csr_matrix(T_ltf_sf), 'constant'),
            EWC_SparseMatrix(sf2.mesh, csr_matrix(T_ltf_tf), 'constant')
        )