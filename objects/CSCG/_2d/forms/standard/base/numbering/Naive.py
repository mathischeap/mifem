# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config.main import *
from screws.freeze.main import FrozenOnly
from tools.linearAlgebra.gathering.regular.chain_matrix.main import Gathering_Matrix, Gathering_Vector



class _2dCSCG_Standard_Form_Numbering_Naive(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()




    def _2dCSCG_0Form_Inner(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._0Form_Inner_no_parameters()
        else:
            raise NotImplementedError()

    def _2dCSCG_0Form_Outer(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._0Form_Outer_no_parameters()
        else:
            raise NotImplementedError()



    def _2dCSCG_1Form_Inner(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._1Form_Inner_no_parameters()
        else:
            raise NotImplementedError()

    def _2dCSCG_1Form_Outer(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._1Form_Outer_no_parameters()
        else:
            raise NotImplementedError()



    def _2dCSCG_2Form_Inner(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._2Form_Inner_no_parameters()
        else:
            raise NotImplementedError()

    def _2dCSCG_2Form_Outer(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._2Form_Outer_no_parameters()
        else:
            raise NotImplementedError()







    def _0Form_Inner_no_parameters(self):
        """
        Do the numbering if it is a standard inner 0-form

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        gathering_matrix = dict()
        element_num = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.num.basis
        extraInfo = None
        if self._sf_.IS.hybrid:
            for i in self._sf_.mesh.elements:
                gathering_matrix[i] = Gathering_Vector(i, range(i * numOfBasis, (i + 1) * numOfBasis))
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if RANK == MASTER_RANK:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1) to make it work for periodic domain."


        if RANK != MASTER_RANK:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            COMM.send([element_map, element_indices], dest=MASTER_RANK, tag=RANK)
            global_numbering = COMM.recv(source=MASTER_RANK, tag=RANK)
        else:
            p = self._sf_.p
            p = (p[0]+1, p[1]+1)
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'DURL' # not an error, this is other side name.
            sidePairDict = {'U':'D', 'D':'U', 'L':'R', 'R':'L'}
            for i in range(SIZE):
                if i == MASTER_RANK:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = COMM.recv(source=i, tag=i)

                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        numberingCache[0] = np.arange(numOfBasis).reshape(p, order='F')
                        currentNumber += numOfBasis
                    else:
                        notNumberedPlaces = numberingCache[k] == -1
                        howManyNotNumbered = len(numberingCache[k][notNumberedPlaces])
                        # assert howManyNotNumbered < numOfBasis, "We make sure at least a part is numbered!"
                        if howManyNotNumbered > 0:
                            numberingCache[k][notNumberedPlaces] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered
                            )
                            currentNumber += howManyNotNumbered

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = -np.ones(p, dtype=int)
                            if otherElement > k:
                                self.___PRIVATE_for_0Form_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )
                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == MASTER_RANK:
                    global_numbering = toBeSentAway
                else:
                    COMM.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = global_numbering[i].ravel('F')
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_0Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):

        if fromSide == 'U'  : data = numberingCache[fromElement][0 , :]
        elif fromSide == 'D': data = numberingCache[fromElement][-1, :]
        elif fromSide == 'L': data = numberingCache[fromElement][ :, 0]
        elif fromSide == 'R': data = numberingCache[fromElement][ :,-1]
        else: raise Exception()

        if toSide == 'U'  : numberingCache[toElement][ 0, :] = data
        elif toSide == 'D': numberingCache[toElement][-1, :] = data
        elif toSide == 'L': numberingCache[toElement][ :, 0] = data
        elif toSide == 'R': numberingCache[toElement][ :,-1] = data
        else: raise Exception()

    def _0Form_Outer_no_parameters(self):
        return self._0Form_Inner_no_parameters()




    def _1Form_Inner_no_parameters(self):
        """
        Do the numbering if it is a standard inner 1-form

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        gathering_matrix = dict()
        element_num = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.num.basis
        extraInfo = None
        if self._sf_.IS.hybrid:
            for i in self._sf_.mesh.elements:
                gathering_matrix[i] = Gathering_Vector(i, range(i * numOfBasis, (i + 1) * numOfBasis))
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if RANK == MASTER_RANK:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1) to make it work for periodic domain."

        if RANK != MASTER_RANK:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            COMM.send([element_map, element_indices], dest=MASTER_RANK, tag=RANK)
            global_numbering = COMM.recv(source=MASTER_RANK, tag=RANK)
        else:
            p = self._sf_.p
            numOfBasisComponents = self._sf_.num.basis_components
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'DURL' # not an error, this is other side name.
            sidePairDict = {'U':'D', 'D':'U', 'L':'R', 'R':'L'}
            for i in range(SIZE):
                if i == MASTER_RANK:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = COMM.recv(source=i, tag=i)
                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        nC0 = np.arange(currentNumber, currentNumber + numOfBasisComponents[0]).reshape(
                            (p[0], p[1] + 1), order='F')
                        currentNumber += numOfBasisComponents[0]
                        nC1 = np.arange(currentNumber, currentNumber + numOfBasisComponents[1]).reshape(
                            (p[0] + 1, p[1]), order='F')
                        currentNumber += numOfBasisComponents[1]
                        numberingCache[0] = [nC0, nC1]
                        del nC0, nC1
                    else:
                        notNumberedPlaces0 = numberingCache[k][0] == -1
                        notNumberedPlaces1 = numberingCache[k][1] == -1
                        howManyNotNumbered0 = len(numberingCache[k][0][notNumberedPlaces0])
                        howManyNotNumbered1 = len(numberingCache[k][1][notNumberedPlaces1])
                        # assert howManyNotNumbered < numOfBasis, "We make sure at least a part is numbered!"
                        if howManyNotNumbered0 > 0:
                            numberingCache[k][0][notNumberedPlaces0] = np.arange(
                                currentNumber, currentNumber + howManyNotNumbered0
                            )
                            currentNumber += howManyNotNumbered0
                        if howManyNotNumbered1 > 0:
                            numberingCache[k][1][notNumberedPlaces1] = np.arange(
                                currentNumber, currentNumber + howManyNotNumbered1
                            )
                            currentNumber += howManyNotNumbered1

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str):  # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = [-np.ones((p[0], p[1] + 1), dtype=int),
                                                                -np.ones((p[0] + 1, p[1]), dtype=int)]
                            if otherElement > k:
                                self.___PRIVATE_for_1FormI_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )
                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == MASTER_RANK:
                    global_numbering = toBeSentAway
                else:
                    COMM.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = np.concatenate([global_numbering[i][0].ravel('F'),
                                     global_numbering[i][1].ravel('F')])
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_1FormI_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):

        data0, data1, data2 = None, None, None

        if fromSide == 'U':
            data1 = numberingCache[fromElement][1][0, :]
        elif fromSide == 'D':
            data1 = numberingCache[fromElement][1][-1, :]
        elif fromSide == 'L':
            data0 = numberingCache[fromElement][0][:, 0]
        elif fromSide == 'R':
            data0 = numberingCache[fromElement][0][:, -1]
        else:
            raise Exception()
        if toSide == 'U':
            numberingCache[toElement][1][0, :] = data1
        elif toSide == 'D':
            numberingCache[toElement][1][-1, :] = data1
        elif toSide == 'L':
            numberingCache[toElement][0][:, 0] = data0
        elif toSide == 'R':
            numberingCache[toElement][0][:, -1] = data0
        else:
            raise Exception()

    def _1Form_Outer_no_parameters(self):
        """
        Do the numbering if it is a standard outer 1-form

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        gathering_matrix = dict()
        element_num = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.num.basis
        extraInfo = None
        if self._sf_.IS.hybrid:
            for i in self._sf_.mesh.elements:
                gathering_matrix[i] = Gathering_Vector(i, range(i * numOfBasis, (i + 1) * numOfBasis))
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if RANK == MASTER_RANK:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1) to make it work for periodic domain."

        if RANK != MASTER_RANK:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            COMM.send([element_map, element_indices], dest=MASTER_RANK, tag=RANK)
            global_numbering = COMM.recv(source=MASTER_RANK, tag=RANK)
        else:
            p = self._sf_.p
            numOfBasisComponents = self._sf_.num.basis_components
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'DURL' # not an error, this is other side name.
            sidePairDict = {'U':'D', 'D':'U', 'L':'R', 'R':'L'}
            for i in range(SIZE):
                if i == MASTER_RANK:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = COMM.recv(source=i, tag=i)
                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        nC0 = np.arange(currentNumber, currentNumber+numOfBasisComponents[0]).reshape(
                            (p[0]+1, p[1]), order='F')
                        currentNumber += numOfBasisComponents[0]
                        nC1 = np.arange(currentNumber, currentNumber+numOfBasisComponents[1]).reshape(
                            (p[0], p[1]+1), order='F')
                        currentNumber += numOfBasisComponents[1]
                        numberingCache[0] = [nC0, nC1]
                        del nC0, nC1
                    else:
                        notNumberedPlaces0 = numberingCache[k][0] == -1
                        notNumberedPlaces1 = numberingCache[k][1] == -1
                        howManyNotNumbered0 = len(numberingCache[k][0][notNumberedPlaces0])
                        howManyNotNumbered1 = len(numberingCache[k][1][notNumberedPlaces1])
                        # assert howManyNotNumbered < numOfBasis, "We make sure at least a part is numbered!"
                        if howManyNotNumbered0 > 0:
                            numberingCache[k][0][notNumberedPlaces0] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered0
                            )
                            currentNumber += howManyNotNumbered0
                        if howManyNotNumbered1 > 0:
                            numberingCache[k][1][notNumberedPlaces1] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered1
                            )
                            currentNumber += howManyNotNumbered1

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = [-np.ones((p[0]+1, p[1]), dtype=int),
                                                                -np.ones((p[0], p[1]+1), dtype=int)]
                            if otherElement > k:
                                self.___PRIVATE_for_1FormO_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )
                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == MASTER_RANK:
                    global_numbering = toBeSentAway
                else:
                    COMM.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = np.concatenate([global_numbering[i][0].ravel('F'),
                                     global_numbering[i][1].ravel('F')])
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_1FormO_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):
        if fromSide == 'U'  : data = numberingCache[fromElement][0][0 , :]
        elif fromSide == 'D': data = numberingCache[fromElement][0][-1, :]
        elif fromSide == 'L': data = numberingCache[fromElement][1][ :, 0]
        elif fromSide == 'R': data = numberingCache[fromElement][1][ :,-1]
        else: raise Exception()

        if toSide == 'U'  : numberingCache[toElement][0][ 0, :] = data
        elif toSide == 'D': numberingCache[toElement][0][-1, :] = data
        elif toSide == 'L': numberingCache[toElement][1][ :, 0] = data
        elif toSide == 'R': numberingCache[toElement][1][ :,-1] = data
        else: raise Exception()




    def _2Form_Inner_no_parameters(self):
        """
        Do the numbering if it is a standard inner 2-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        gathering_matrix = dict()
        element_num = self._sf_.mesh.elements.num
        numOfBasis = self._sf_.num.basis
        extraInfo = None
        for i in self._sf_.mesh.elements:
            gathering_matrix[i] = Gathering_Vector(i, range(i * numOfBasis, (i + 1) * numOfBasis))
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_2dCSCG')
        numOfDofs = numOfBasis * element_num
        return gathering_matrix, numOfDofs, extraInfo

    def _2Form_Outer_no_parameters(self):
        return self._2Form_Inner_no_parameters()