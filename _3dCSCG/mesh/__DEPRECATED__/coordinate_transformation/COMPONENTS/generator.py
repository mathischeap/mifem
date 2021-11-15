# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Fri May  3 09:47:58 2019
Aerodynamics, AE
TU Delft
"""
from SCREWS.frozen import FrozenOnly

class CoordinateTransformationGenerators(FrozenOnly):
    """ Here we group all generators for class CoordinateTransformation3D."""
    def __init__(self, ct3):
        """ """
        self._ct3_ = ct3
        self._freeze_self_()
    
    def mapping(self, evaluation_points):
        return self._ct3_.___mapping_generator___(evaluation_points)
    
    def Jacobian_matrix(self, evaluation_points):
        return self._ct3_.___Jacobian_matrix_generator___(evaluation_points)
        
    def Jacobian(self, evaluation_points):
        return self._ct3_.___Jacobian_generator___(evaluation_points)
    
    def inverse_Jacobian_matrix(self, evaluation_points):
        return self._ct3_.___inverse_Jacobian_matrix_generator___(evaluation_points)
    
    def inverse_metric_matrix(self, evaluation_points):
        return self._ct3_.___inverse_metric_matrix_generator___(evaluation_points)