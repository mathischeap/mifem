# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Thu May 23 11:07:25 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
import numpy as np
from screws.frozen import FrozenOnly

class CTEXTTBase(FrozenOnly):
    """ 
    Parent of CoordinateTransformationTrace3D. 
    
    In Trace, the very import difference is that we take NOT mesh-grid inputs
    to evaluate the mapping and so on. This is because we have different faces
    (3D) or edges (2D) for the trace mapping. If we use the mesh-grid points
    as we did in CoordinateTransformation, then, for example, in 3D case, we 
    needs 6 inputs of shape = (2,), at least, which is not a smart way.
    
    """
    def __init__(self, ct):
        """ """
        self._ct_ = ct
#        self._mapping_ = None
#        self._Jacobian_matrix_ = None
#        self._metric_matrix_ = None
        self._freeze_self_()
        
#    def _reset_(self):
#        self._mapping_ = None
#        self._Jacobian_matrix_ = None
#        self._metric_matrix_ = None
    
    @property
    def _mesh_(self):
        return self._ct_._mesh_
    
    @property
    def ndim(self):
        """ 
        this is a trace mapping is in n dimensonal object; itself is a n-1 dimensional
        one.
        
        """
        return self._ct_.ndim
    
    def ___generate_trace_evaluation_points___(self):
        """ 
        When even we try to compute the trace mapping or trace Jacobian_matrix,
        we run this method before hands to generate proper points (in reference
        coordinates) for 4 edges(2D) or 6 sides(3D).
        
        This looks very bad, since if we have done trace mapping, when we 
        further compute trace Jacobian, we will repeat it again, why not just 
        store it? No, we do not do this, because we always search 
        evaluation_points_gird from `self._ct_`, and when we reset 
        evaluation_points_gird, we will have to reset the stored value as well.
        Of course, this is doable, but I think it makes the code a little bit
        more un-readable. And what is more, this process is very fast, so, who
        cares if we run it one more time?
        
        """
        raise Exception(" <CoordinateTransformation.Trace> : To be overwritten.")
    
    @property
    def mapping(self):
        """ 
        The mapping. To compute it we just need to employ the 
        CoordinateTransformation.
        
        Returns
        -------
        self._mapping_ : dict
            Unlike the CoordinateTransformation.mapping which must be of 
            structured data sturcture (so we put it in a ndarray), we here put
            it in a dict just like what we have in meshComponents.trace.
        
        """
        raise Exception(" <CoordinateTransformation.Trace> : To be overwritten.")
        
    @property
    def Jacobian_matrix(self):
        """ 
        The Jacobian matrix. To compute it we just need to employ the 
        CoordinateTransformation.
        
        Returns
        -------
        self._Jacobian_matrix_ : dict
            As self.mapping, here we also put it in a dict whose keys represent
            the numbering of the trace element. Just like what we always have 
            in meshComponents.trace.
            
        """
        raise Exception(" <CoordinateTransformation.Trace> : To be overwritten.")
    
    @property
    def metric_matrix(self):
        """ The entries of metric_matrix is normally denoted as g_{i,j}. """
#        if self._metric_matrix_ is None:
#            J = self.Jacobian_matrix
#            G = {}
#            for k in self._mesh_.trace.elements.position_representive:
#                Gk = [[None for i in range(self.ndim)] for j in range(self.ndim)]
#                for i in range(self.ndim):
#                    for j in range(i, self.ndim):
#                        Gk[i][j] = J[k][0][i] * J[k][0][j]
#                        for l in range(1, self._ct_.ndim):
#                            Gk[i][j] += J[k][l][i] * J[k][l][j]
#                        if i != j:
#                            Gk[j][i] = Gk[i][j]
#                G[k] = np.array(Gk)
#            self._metric_matrix_ = G
#        return self._metric_matrix_
        J = self.Jacobian_matrix
        G = {}
        for k in self._mesh_.trace.elements.position_representive:
            Gk = [[None for _ in range(self.ndim-1)] for _ in range(self.ndim-1)]
            for i in range(self.ndim-1):
                for j in range(i, self.ndim-1):
                    Gk[i][j] = J[k][0][i] * J[k][0][j]
                    for l in range(1, self.ndim):
                        Gk[i][j] += J[k][l][i] * J[k][l][j]
                    if i != j:
                        Gk[j][i] = Gk[i][j]
            G[k] = np.array(Gk)
        return G
