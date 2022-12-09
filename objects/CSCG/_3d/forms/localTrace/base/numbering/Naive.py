# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:30 PM
"""
from components.freeze.main import FrozenOnly

from tools.elementwiseCache.gathering.vector import Gathering_Vector
from tools.elementwiseCache.gathering.regular.matrix.main import Gathering_Matrix




class _3dCSCG_LocalTrace_Form_Numbering_Naive(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def _3dCSCG_0LocalTrace(self):
        """Do the numbering if it is a local-trace-0-form:
        :class:`_3dCSCG.form.localTrace._0ltf._0LocalTrace`.

        :returns: A tuple of 3 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...) -- Extra numbering information.
        """
        if self._ltf_.numbering._parameters_ == dict():
            return self._no_parameter_numbering()
        else:
            raise NotImplementedError()

    def _no_parameter_numbering(self):
        """"""
        mesh = self._ltf_.mesh
        gathering_matrix = dict()
        element_num = mesh.elements.num
        numOfBasis = self._ltf_.num.basis

        extraInfo = None
        for i in self._ltf_.mesh.elements:
            gathering_matrix[i] = Gathering_Vector(i, range(i * numOfBasis, (i + 1) * numOfBasis))
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
        numOfDofs = numOfBasis * element_num
        return gathering_matrix, numOfDofs, extraInfo

    def _3dCSCG_1LocalTrace(self):
        """Do the numbering if it is a local-trace-1-form:
        :class:`_3dCSCG.form.localTrace._1ltf._1LocalTrace`.

        :returns: A tuple of 3 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...) -- Extra numbering information.
        """
        if self._ltf_.numbering._parameters_ == dict():
            return self._no_parameter_numbering()
        else:
            raise NotImplementedError()

    def _3dCSCG_2LocalTrace(self):
        """Do the numbering if it is a local-trace-2-form:
        :class:`_3dCSCG.form.localTrace._2ltf._2LocalTrace`.

        :returns: A tuple of 3 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...) -- Extra numbering information.
        """
        if self._ltf_.numbering._parameters_ == dict():
            return self._no_parameter_numbering()
        else:
            raise NotImplementedError()