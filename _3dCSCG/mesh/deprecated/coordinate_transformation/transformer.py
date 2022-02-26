# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Wed Sep 12 14:39:12 2018
Aerodynamics, AE
TU Delft
"""
from _3dCSCG.mesh.deprecated.coordinate_transformation.structuredMeshCTBase.base import CTBase
from _3dCSCG.mesh.deprecated.coordinate_transformation.COMPONENTS.dump import CoordinateTransformationDump
from _3dCSCG.mesh.deprecated.coordinate_transformation.COMPONENTS.trace import CoordinateTransformationTrace
# from _3dCSCG.mesh.__DEPRECATED__.coordinate_transformation.COMPONENTS.generator import CoordinateTransformationGenerators

class CoordinateTransformation(CTBase):
    """ The class CoordinateTransformation for 3D. """
    def __init__(self, mesh):
        self._trace_ = CoordinateTransformationTrace(self)
        self._dump_ = CoordinateTransformationDump(self)
        # self._generators_ = CoordinateTransformationGenerators(self)
        super().__init__(mesh)
    
    @property
    def trace(self):
        """ """
        return self._trace_
    
    @property
    def dump(self):
        """ """
        return self._dump_
    
    # @property
    # def generators(self):
    #     return self._generators_
    
    @property
    def Jacobian(self):
        """ Determinant of the Jacobian matrix. """
        J = self.Jacobian_matrix
        return + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]\
               - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]
    
    def ___Jacobian_generator___(self, evaluation_points):
        """ """
        i = 0
        while i < self._mesh_.elements.num:
            J = self.___compute_Jacobian_matrix_i___(evaluation_points, i)
            yield + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]\
                  - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]
            i += 1
    
    @property
    def inverse_Jacobian(self):
        """ Determinant of the inverse Jacobian matrix. """
        iJ = self.inverse_Jacobian_matrix
        return + iJ[0][0]*iJ[1][1]*iJ[2][2] + iJ[0][1]*iJ[1][2]*iJ[2][0] + iJ[0][2]*iJ[1][0]*iJ[2][1]\
               - iJ[0][0]*iJ[1][2]*iJ[2][1] - iJ[0][1]*iJ[1][0]*iJ[2][2] - iJ[0][2]*iJ[1][1]*iJ[2][0]
        
    @property
    def inverse_Jacobian_matrix(self):
        """ The inverse Jacobian matrix. """
        J = self.Jacobian_matrix
        Jacobian = + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]\
                   - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]
        reciprocal_Jacobian = 1 / Jacobian
        del Jacobian
        iJ00 = reciprocal_Jacobian * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
        iJ01 = reciprocal_Jacobian * (J[2][1]*J[0][2] - J[2][2]*J[0][1])
        iJ02 = reciprocal_Jacobian * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
        iJ10 = reciprocal_Jacobian * (J[1][2]*J[2][0] - J[1][0]*J[2][2])
        iJ11 = reciprocal_Jacobian * (J[2][2]*J[0][0] - J[2][0]*J[0][2])
        iJ12 = reciprocal_Jacobian * (J[0][2]*J[1][0] - J[0][0]*J[1][2])
        iJ20 = reciprocal_Jacobian * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
        iJ21 = reciprocal_Jacobian * (J[2][0]*J[0][1] - J[2][1]*J[0][0])
        iJ22 = reciprocal_Jacobian * (J[0][0]*J[1][1] - J[0][1]*J[1][0])
        return [[iJ00, iJ01, iJ02],
                [iJ10, iJ11, iJ12],
                [iJ20, iJ21, iJ22]]
    
    def ___inverse_Jacobian_matrix_generator___(self, evaluation_points):
        """ """
        detJg = self.___Jacobian_generator___(evaluation_points)
        i = 0
        while i < self._mesh_.elements.num:
            detJ = 1 / detJg.__next__()
            J = self.___compute_Jacobian_matrix_i___(evaluation_points, i)
            iJ00 = detJ * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
            iJ01 = detJ * (J[2][1]*J[0][2] - J[2][2]*J[0][1])
            iJ02 = detJ * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
            iJ10 = detJ * (J[1][2]*J[2][0] - J[1][0]*J[2][2])
            iJ11 = detJ * (J[2][2]*J[0][0] - J[2][0]*J[0][2])
            iJ12 = detJ * (J[0][2]*J[1][0] - J[0][0]*J[1][2])
            iJ20 = detJ * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
            iJ21 = detJ * (J[2][0]*J[0][1] - J[2][1]*J[0][0])
            iJ22 = detJ * (J[0][0]*J[1][1] - J[0][1]*J[1][0])
            yield ((iJ00, iJ01, iJ02),
                   (iJ10, iJ11, iJ12),
                   (iJ20, iJ21, iJ22))
            i += 1

    @staticmethod
    def ___method___(J):
        """
        If we already knwo J, we can use this method to compute the rest without 
        computing J again. This saves a bit time.
        
        Parameters
        ----------
        J :
            Jacobian matrix
            
        Returns
        -------
        output1: 
            The inverse_Jacobian_matrix.
        Jacobian :
            Jacobian.
        reciprocal_Jacobian :
            1 / Jacobian
        
        """
        Jacobian = + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]\
                   - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]
        reciprocal_Jacobian = 1 / Jacobian
        iJ00 = reciprocal_Jacobian * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
        iJ01 = reciprocal_Jacobian * (J[2][1]*J[0][2] - J[2][2]*J[0][1])
        iJ02 = reciprocal_Jacobian * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
        iJ10 = reciprocal_Jacobian * (J[1][2]*J[2][0] - J[1][0]*J[2][2])
        iJ11 = reciprocal_Jacobian * (J[2][2]*J[0][0] - J[2][0]*J[0][2])
        iJ12 = reciprocal_Jacobian * (J[0][2]*J[1][0] - J[0][0]*J[1][2])
        iJ20 = reciprocal_Jacobian * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
        iJ21 = reciprocal_Jacobian * (J[2][0]*J[0][1] - J[2][1]*J[0][0])
        iJ22 = reciprocal_Jacobian * (J[0][0]*J[1][1] - J[0][1]*J[1][0])
        return ((iJ00, iJ01, iJ02),
                (iJ10, iJ11, iJ12),
                (iJ20, iJ21, iJ22)), Jacobian, reciprocal_Jacobian