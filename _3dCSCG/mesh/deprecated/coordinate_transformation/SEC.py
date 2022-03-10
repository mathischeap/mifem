# -*- coding: utf-8 -*-
"""
It is not like `.trace`, `.dump` or `.generator` that are attributes of the `ct`, here
we actually wrap a `ct` and have a new object.

Yi Zhang (C)
Created on Fri May  3 09:47:29 2019
Aerodynamics, AE
TU Delft
"""
from screws.freeze.main import FrozenOnly


class CoordinateTransformationSingleElementComputer(FrozenOnly):
    """ Here we store the methods those can be called for single element."""
    def __init__(self, ct3):
        self._ct3_ = ct3
        assert self._ct3_.evaluation_points is not None, \
            " <CoordinateTransformation3DSingleElementComputer> : first put evaluation_points for the CT3 object."
        self._evl_pts_ = self._ct3_.evaluation_points
        self._freeze_self_()
        
    @property
    def evaluation_points(self):
        return self._evl_pts_
    
    def mapping(self, i):
        """ Compute the mapping for `i`th element."""
        return self._ct3_.___compute_mapping_element_i___(self.evaluation_points, i)
    
    def Jacobian_matrix(self, i):
        """ Compute the Jacobian_matrix for `i`th element."""
        return self._ct3_.___compute_Jacobian_matrix_i___(self.evaluation_points, i)
    
    def Jacobian(self, i):
        J = self.Jacobian_matrix(i)
        return + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]\
               - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]
        
    def inverse_Jacobian_matrix(self, i):
        detJ = 1 / self.Jacobian(i)
        J = self.Jacobian_matrix(i)
        iJ00 = detJ * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
        iJ01 = detJ * (J[2][1]*J[0][2] - J[2][2]*J[0][1])
        iJ02 = detJ * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
        iJ10 = detJ * (J[1][2]*J[2][0] - J[1][0]*J[2][2])
        iJ11 = detJ * (J[2][2]*J[0][0] - J[2][0]*J[0][2])
        iJ12 = detJ * (J[0][2]*J[1][0] - J[0][0]*J[1][2])
        iJ20 = detJ * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
        iJ21 = detJ * (J[2][0]*J[0][1] - J[2][1]*J[0][0])
        iJ22 = detJ * (J[0][0]*J[1][1] - J[0][1]*J[1][0])
        return ((iJ00, iJ01, iJ02),
                (iJ10, iJ11, iJ12),
                (iJ20, iJ21, iJ22))
    
    def ___inverse_Jacobian_matrix_ep___(self, i, ep):
        J = self._ct3_.___compute_Jacobian_matrix_i___(ep, i)
        detJ = 1/(+ J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] + J[0][2]*J[1][0]*J[2][1]
                  - J[0][0]*J[1][2]*J[2][1] - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0])
        iJ00 = detJ * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
        iJ01 = detJ * (J[2][1]*J[0][2] - J[2][2]*J[0][1])
        iJ02 = detJ * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
        iJ10 = detJ * (J[1][2]*J[2][0] - J[1][0]*J[2][2])
        iJ11 = detJ * (J[2][2]*J[0][0] - J[2][0]*J[0][2])
        iJ12 = detJ * (J[0][2]*J[1][0] - J[0][0]*J[1][2])
        iJ20 = detJ * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
        iJ21 = detJ * (J[2][0]*J[0][1] - J[2][1]*J[0][0])
        iJ22 = detJ * (J[0][0]*J[1][1] - J[0][1]*J[1][0])
        return ((iJ00, iJ01, iJ02),
                (iJ10, iJ11, iJ12),
                (iJ20, iJ21, iJ22))