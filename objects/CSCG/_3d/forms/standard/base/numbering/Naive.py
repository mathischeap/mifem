# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config.main import *
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.gathering.chain_matrix.main import Gathering_Matrix, Gathering_Vector



class _3dCSCG_Standard_Form_Numbering_Naive(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._mesh_ = sf.mesh
        self._freeze_self_()





    def _3dCSCG_0Form(self):
        """
        Do the numbering if it is a standard 0-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        parameters = self._sf_.numbering._parameters_
        if self._sf_.numbering._parameters_ == dict():
            return self._0Form_no_parameters()
        elif len(parameters) == 1 and 'SingleRegionSideCrack' in parameters:
            # have one regions side crack ...
            assert not self._sf_.IS.hybrid, \
                "0-form can not be hybrid for single-regions-side-crack-numbering."
            return self._0Form_a_region_side_crack()
        else:
            raise NotImplementedError()


    def _3dCSCG_1Form(self):
        """
        Do the numbering if it is a standard 1-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._1Form_no_parameters()
        else:
            raise NotImplementedError()


    def _3dCSCG_2Form(self):
        """
        Do the numbering if it is a standard 2-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._2Form_no_parameters()
        else:
            raise NotImplementedError()


    def _3dCSCG_3Form(self):
        """
        Do the numbering if it is a standard 3-form.

        :returns: A tuple of 3 outputs:

            1. (str) -- The global numbering in local elements.
            2. (int) -- Number of dofs in this core.
            3. (None,...)-- Extra numbering information.
        """
        if self._sf_.numbering._parameters_ == dict():
            return self._3Form_no_parameters()
        else:
            raise NotImplementedError()






    def _3Form_no_parameters(self):
        """
        Do the numbering if it is a standard 3-form.

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
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
        numOfDofs = numOfBasis * element_num
        return gathering_matrix, numOfDofs, extraInfo





    def _2Form_no_parameters(self):
        """
        Do the numbering if it is a standard 2-form:
        :class:`_3dCSCG.form.standard._2_form._2Form`.

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
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...


        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if rAnk == mAster_rank:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1, >1) to make it work for periodic domain."

        if rAnk != mAster_rank:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            cOmm.send([element_map, element_indices], dest=mAster_rank, tag=rAnk)
            global_numbering = cOmm.recv(source=mAster_rank, tag=rAnk)
        else:
            p = self._sf_.p
            numOfBasisComponents = self._sf_.num.basis_components
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'SNEWFB' # not an error, this is other side name.
            sidePairDict = {'N':'S', 'S':'N', 'W':'E', 'E':'W', 'B':'F', 'F':'B'}
            for i in range(sIze):
                if i == mAster_rank:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = cOmm.recv(source=i, tag=i)
                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        nC0 = np.arange(currentNumber, currentNumber+numOfBasisComponents[0]).reshape(
                            (p[0]+1, p[1], p[2]), order='F')
                        currentNumber += numOfBasisComponents[0]
                        nC1 = np.arange(currentNumber, currentNumber+numOfBasisComponents[1]).reshape(
                            (p[0], p[1]+1, p[2]), order='F')
                        currentNumber += numOfBasisComponents[1]
                        nC2 = np.arange(currentNumber, currentNumber+numOfBasisComponents[2]).reshape(
                            (p[0], p[1], p[2]+1), order='F')
                        currentNumber += numOfBasisComponents[2]
                        numberingCache[0] = [nC0, nC1, nC2]
                        del nC0, nC1, nC2
                    else:
                        notNumberedPlaces0 = numberingCache[k][0] == -1
                        notNumberedPlaces1 = numberingCache[k][1] == -1
                        notNumberedPlaces2 = numberingCache[k][2] == -1
                        howManyNotNumbered0 = len(numberingCache[k][0][notNumberedPlaces0])
                        howManyNotNumbered1 = len(numberingCache[k][1][notNumberedPlaces1])
                        howManyNotNumbered2 = len(numberingCache[k][2][notNumberedPlaces2])
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
                        if howManyNotNumbered2 > 0:
                            numberingCache[k][2][notNumberedPlaces2] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered2
                            )
                            currentNumber += howManyNotNumbered2

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = [-np.ones((p[0]+1, p[1], p[2]), dtype=int),
                                                                -np.ones((p[0], p[1]+1, p[2]), dtype=int),
                                                                -np.ones((p[0], p[1], p[2]+1), dtype=int)]
                            if otherElement > k:
                                self.___PRIVATE_for_2Form_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )
                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == mAster_rank:
                    global_numbering = toBeSentAway
                else:
                    cOmm.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = np.concatenate([global_numbering[i][0].ravel('F'),
                                     global_numbering[i][1].ravel('F'),
                                     global_numbering[i][2].ravel('F')])
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_2Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):
        if fromSide == 'N'  : data = numberingCache[fromElement][0][0 , :, :]
        elif fromSide == 'S': data = numberingCache[fromElement][0][-1, :, :]
        elif fromSide == 'W': data = numberingCache[fromElement][1][ :, 0, :]
        elif fromSide == 'E': data = numberingCache[fromElement][1][ :,-1, :]
        elif fromSide == 'B': data = numberingCache[fromElement][2][ :, :, 0]
        elif fromSide == 'F': data = numberingCache[fromElement][2][ :, :,-1]
        else: raise Exception()

        if toSide == 'N'  : numberingCache[toElement][0][ 0, :, :] = data
        elif toSide == 'S': numberingCache[toElement][0][-1, :, :] = data
        elif toSide == 'W': numberingCache[toElement][1][ :, 0, :] = data
        elif toSide == 'E': numberingCache[toElement][1][ :,-1, :] = data
        elif toSide == 'B': numberingCache[toElement][2][ :, :, 0] = data
        elif toSide == 'F': numberingCache[toElement][2][ :, :,-1] = data
        else: raise Exception()




    def _1Form_no_parameters(self):
        """
        Do the numbering if it is a standard 1-form:
        :class:`_3dCSCG.form.standard._1_form._1Form`.

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
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if rAnk == mAster_rank:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1, >1) to make it work for periodic domain."

        if rAnk != mAster_rank:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            cOmm.send([element_map, element_indices], dest=mAster_rank, tag=rAnk)
            global_numbering = cOmm.recv(source=mAster_rank, tag=rAnk)
        else:
            p = self._sf_.p
            numOfBasisComponents = self._sf_.num.basis_components
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'SNEWFB' # not an error, this is other side name.
            sidePairDict = {'N':'S', 'S':'N', 'W':'E', 'E':'W', 'B':'F', 'F':'B'}
            for i in range(sIze):
                if i == mAster_rank:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = cOmm.recv(source=i, tag=i)
                for k in element_indices:
                    if k == 0:
                        assert numberingCache == dict(), "We make sure #0 element is numbered firstly."
                        nC0 = np.arange(currentNumber, currentNumber+numOfBasisComponents[0]).reshape(
                            (p[0], p[1]+1, p[2]+1), order='F')
                        currentNumber += numOfBasisComponents[0]
                        nC1 = np.arange(currentNumber, currentNumber+numOfBasisComponents[1]).reshape(
                            (p[0]+1, p[1], p[2]+1), order='F')
                        currentNumber += numOfBasisComponents[1]
                        nC2 = np.arange(currentNumber, currentNumber+numOfBasisComponents[2]).reshape(
                            (p[0]+1, p[1]+1, p[2]), order='F')
                        currentNumber += numOfBasisComponents[2]
                        numberingCache[0] = [nC0, nC1, nC2]
                        del nC0, nC1, nC2
                    else:
                        notNumberedPlaces0 = numberingCache[k][0] == -1
                        notNumberedPlaces1 = numberingCache[k][1] == -1
                        notNumberedPlaces2 = numberingCache[k][2] == -1
                        howManyNotNumbered0 = len(numberingCache[k][0][notNumberedPlaces0])
                        howManyNotNumbered1 = len(numberingCache[k][1][notNumberedPlaces1])
                        howManyNotNumbered2 = len(numberingCache[k][2][notNumberedPlaces2])
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
                        if howManyNotNumbered2 > 0:
                            numberingCache[k][2][notNumberedPlaces2] = np.arange(
                                currentNumber, currentNumber+howManyNotNumbered2
                            )
                            currentNumber += howManyNotNumbered2

                    for j, EMki in enumerate(element_map[k]):
                        if isinstance(EMki, str): # on domain boundary
                            pass
                        else:
                            otherElement = EMki
                            otherSide = other_side_name[j]
                            if otherElement not in numberingCache and otherElement > k:
                                numberingCache[otherElement] = [-np.ones((p[0], p[1]+1, p[2]+1), dtype=int),
                                                                -np.ones((p[0]+1, p[1], p[2]+1), dtype=int),
                                                                -np.ones((p[0]+1, p[1]+1, p[2]), dtype=int)]
                            if otherElement > k:
                                self.___PRIVATE_for_1Form_pass_element_side_numbering_from_to___(
                                    k, sidePairDict[otherSide], otherElement, otherSide, numberingCache
                                )
                toBeSentAway = dict()
                for k in element_indices: toBeSentAway[k] = numberingCache[k]
                if i == mAster_rank:
                    global_numbering = toBeSentAway
                else:
                    cOmm.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = np.concatenate([global_numbering[i][0].ravel('F'),
                                     global_numbering[i][1].ravel('F'),
                                     global_numbering[i][2].ravel('F')])
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_1Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):

        data0, data1, data2 = None, None, None

        if fromSide == 'N'  :
            data1 = numberingCache[fromElement][1][0 , :, :]
            data2 = numberingCache[fromElement][2][0 , :, :]
        elif fromSide == 'S':
            data1 = numberingCache[fromElement][1][-1, :, :]
            data2 = numberingCache[fromElement][2][-1, :, :]
        elif fromSide == 'W':
            data0 = numberingCache[fromElement][0][ :, 0, :]
            data2 = numberingCache[fromElement][2][ :, 0, :]
        elif fromSide == 'E':
            data0 = numberingCache[fromElement][0][ :,-1, :]
            data2 = numberingCache[fromElement][2][ :,-1, :]
        elif fromSide == 'B':
            data0 = numberingCache[fromElement][0][ :, :, 0]
            data1 = numberingCache[fromElement][1][ :, :, 0]
        elif fromSide == 'F':
            data0 = numberingCache[fromElement][0][ :, :,-1]
            data1 = numberingCache[fromElement][1][ :, :,-1]
        else: raise Exception()
        if toSide == 'N'  :
            numberingCache[toElement][1][ 0, :, :] = data1
            numberingCache[toElement][2][ 0, :, :] = data2
        elif toSide == 'S':
            numberingCache[toElement][1][-1, :, :] = data1
            numberingCache[toElement][2][-1, :, :] = data2
        elif toSide == 'W':
            numberingCache[toElement][0][ :, 0, :] = data0
            numberingCache[toElement][2][ :, 0, :] = data2
        elif toSide == 'E':
            numberingCache[toElement][0][ :,-1, :] = data0
            numberingCache[toElement][2][ :,-1, :] = data2
        elif toSide == 'B':
            numberingCache[toElement][0][ :, :, 0] = data0
            numberingCache[toElement][1][ :, :, 0] = data1
        elif toSide == 'F':
            numberingCache[toElement][0][ :, :,-1] = data0
            numberingCache[toElement][1][ :, :,-1] = data1
        else: raise Exception()




    def _0Form_no_parameters(self):
        """
        Do the numbering if it is a standard 0-form:
        :class:`_3dCSCG.form.standard._0_form._0Form`.

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
            gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
            numOfDofs = numOfBasis * element_num
            return gathering_matrix, numOfDofs, extraInfo

        global_numbering = None
        # non-hybrid numbering ...
        mesh = self._sf_.mesh
        if mesh.domain.IS.periodic:
            if rAnk == mAster_rank:
                baseElementLayout = mesh.elements.layout
                for rn in baseElementLayout:
                    region = mesh.domain.regions[rn]
                    if region.IS.periodic_to_self:
                        regionElementLayout = baseElementLayout[rn]
                        assert all(np.array(regionElementLayout) > 1), \
                            f" elements.layout[{rn}]={regionElementLayout} wrong," \
                            f" needs (>1, >1, >1) to make it work for periodic domain."

        if rAnk != mAster_rank:
            element_map = mesh.elements.map
            element_indices = mesh.elements.indices
            cOmm.send([element_map, element_indices], dest=mAster_rank, tag=rAnk)
            global_numbering = cOmm.recv(source=mAster_rank, tag=rAnk)
        else:
            p = self._sf_.p
            p = (p[0]+1, p[1]+1, p[2]+1)
            numberingCache = dict()
            currentNumber = 0
            other_side_name = 'SNEWFB' # not an error, this is other side name.
            sidePairDict = {'N':'S', 'S':'N', 'W':'E', 'E':'W', 'B':'F', 'F':'B'}
            for i in range(sIze):
                if i == mAster_rank:
                    element_map = mesh.elements.map
                    element_indices = mesh.elements.indices
                else:
                    element_map, element_indices = cOmm.recv(source=i, tag=i)

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
                if i == mAster_rank:
                    global_numbering = toBeSentAway
                else:
                    cOmm.send(toBeSentAway, dest=i, tag=i)
                for k in element_indices: del numberingCache[k]

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            vector = global_numbering[i].ravel('F')
            gathering_matrix[i] = Gathering_Vector(i, vector)
            dofsPOOL.update(vector)
        gathering_matrix = Gathering_Matrix(gathering_matrix, mesh_type='_3dCSCG')
        numOfDofs = len(dofsPOOL)

        return gathering_matrix, numOfDofs, extraInfo

    @staticmethod
    def ___PRIVATE_for_0Form_pass_element_side_numbering_from_to___(
        fromElement, fromSide, toElement, toSide, numberingCache):
        if fromSide == 'N'  : data = numberingCache[fromElement][0 , :, :]
        elif fromSide == 'S': data = numberingCache[fromElement][-1, :, :]
        elif fromSide == 'W': data = numberingCache[fromElement][ :, 0, :]
        elif fromSide == 'E': data = numberingCache[fromElement][ :,-1, :]
        elif fromSide == 'B': data = numberingCache[fromElement][ :, :, 0]
        elif fromSide == 'F': data = numberingCache[fromElement][ :, :,-1]
        else: raise Exception()

        if toSide == 'N'  : numberingCache[toElement][ 0, :, :] = data
        elif toSide == 'S': numberingCache[toElement][-1, :, :] = data
        elif toSide == 'W': numberingCache[toElement][ :, 0, :] = data
        elif toSide == 'E': numberingCache[toElement][ :,-1, :] = data
        elif toSide == 'B': numberingCache[toElement][ :, :, 0] = data
        elif toSide == 'F': numberingCache[toElement][ :, :,-1] = data
        else: raise Exception()

    def _0Form_a_region_side_crack(self):
        """Do the numbering for 0Form with one regions side crack."""
        crack = self._sf_.numbering._parameters_['SingleRegionSideCrack']
        assert isinstance(crack, str), "parameter of a crack must be a string."
        gathering_matrix, numOfDofs, extraInfo = self._0Form_no_parameters()

        side1, side2 = crack.split('|')
        r1, s1 = side1.split('-')
        s2, r2 = side2.split('-')
        assert r1 in self._mesh_.domain.regions and r2 in self._mesh_.domain.regions
        if s1 == 'N': assert s2 == 'S'
        elif s1 == 'S': assert s2 == 'N'
        elif s1 == 'W': assert s2 == 'E'
        elif s1 == 'E': assert s2 == 'W'
        elif s1 == 'B': assert s2 == 'F'
        elif s1 == 'F': assert s2 == 'B'
        else:
            raise Exception()

        DOFs_2b_updated = set()

        elements_1 = self._mesh_.do.FIND_element_attach_to_region_side(r1, s1)
        elements_2 = self._mesh_.do.FIND_element_attach_to_region_side(r2, s2)

        assert elements_1.shape == elements_2.shape
        I, J = elements_1.shape
        involved_elements = list()
        for j in range(J):
            for i in range(I):
                dofs_1, dofs_2 = None, None

                if saFe_mode:
                    e1 = elements_1[i, j]
                    if e1 in self._mesh_.elements:
                        dofs_1 = self._sf_.numbering.___PRIVATE_DO_find_dofs_on_element_side___(e1, s1, GM=gathering_matrix)
                        # dofs_1 = self._sf_.numbering.do.find.dofs_on_element_side(e1, s1, GM=gathering_matrix)

                e2 = elements_2[i, j]
                if e2 in self._mesh_.elements:
                    involved_elements.append(e2)
                    dofs_2 = self._sf_.numbering.___PRIVATE_DO_find_dofs_on_element_side___(e2, s2, GM=gathering_matrix)

                if saFe_mode:
                    if dofs_2 is not None and dofs_1 is not None:
                        np.testing.assert_array_equal(dofs_1, dofs_2)

                if dofs_2 is not None:
                    DOFs_2b_updated.update(dofs_2)
                else:
                    pass

        GLOBAL_num_dofs = gathering_matrix.GLOBAL_num_dofs

        ALL_DOFs_2b_updated = cOmm.gather(DOFs_2b_updated, root=mAster_rank)
        if rAnk == mAster_rank:
            DOFS = set()
            for _ in ALL_DOFs_2b_updated:
                DOFS.update(_)
            S = sorted(DOFS)
            all_new_number_dict = dict()
            for i, old_number in enumerate(S):
                all_new_number_dict[old_number] = GLOBAL_num_dofs + i
        else:
            all_new_number_dict = None
        all_new_number_dict = cOmm.bcast(all_new_number_dict, root=mAster_rank)
        PAIRS = dict()
        for i in DOFs_2b_updated:
            PAIRS[i] = all_new_number_dict[i]
        del all_new_number_dict
        extraInfo = {f'DOF_PAIRS of SingleRegionSideCrack': PAIRS}

        GVD = dict()
        for i in self._mesh_.elements:
            if i in involved_elements:
                old_fv = gathering_matrix[i].full_vector
                indices = self._sf_.numbering._localSideCache0_[s2]
                for ind in indices: old_fv[ind] = PAIRS[old_fv[ind]]
                GVD[i] = Gathering_Vector(i, old_fv)
            else:
                GVD[i] = gathering_matrix[i]

        del gathering_matrix
        new_gathering_matrix = Gathering_Matrix(GVD, mesh_type='_3dCSCG')

        dofsPOOL = set()
        for i in self._sf_.mesh.elements:
            dofsPOOL.update(new_gathering_matrix[i].full_vector)
        numOfDofs = len(dofsPOOL)

        return new_gathering_matrix, numOfDofs, extraInfo














if __name__ == '__main__':
    pass