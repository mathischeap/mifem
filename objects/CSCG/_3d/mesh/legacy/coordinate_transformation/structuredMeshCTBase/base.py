# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Wed May 22 20:53:12 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
import numpy as np
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.legacy.coordinate_transformation.structuredMeshCTBase.MODULES.methods import CTMODMethods
from objects.CSCG._3d.mesh.legacy.coordinate_transformation.structuredMeshCTBase.MODULES.generators import CTMODGenerators

class CTBase(CTMODGenerators, CTMODMethods, FrozenOnly):
    """ The parent class CoordinateTransformation. """
    def __init__(self, mesh):
        """ """
        self._mesh_ = mesh
#        self._mapping_ = None
#        self._Jacobian_matrix_ = None
        self._evaluation_points_ = None
        self._evaluation_points_grid_ = None
        self._freeze_self_()
    
    def __call__(self):
        return self.mapping
    
    @property
    def ndim(self):
        return self._mesh_.ndim
    
    def _reset_(self):
#        self.trace._reset_()
#        self._mapping_ = None
#        self._Jacobian_matrix_ = None
        self._evaluation_points_ = None
        self._evaluation_points_grid_ = None
        
    #________________ evaluation points________________________________________________
    @property
    def evaluation_points(self):
        return self._evaluation_points_
    
    @property
    def evaluation_points_grid(self):
        return self._evaluation_points_grid_
    
    def evaluated_at(self, *args):
        self._reset_()
        self._set_evaluation_points_(*args)
    
    def evaluated_at_meshgrid(self, *args):
        """ 
        When we want to evaluate the mapping by giving 1-D data, we will do the
        'meshgrid' (no ravel will do to it further) to it later.
        
        """
        self._reset_()
        self._set_evaluation_points_grid_(*args)
    
    def _set_evaluation_points_grid_(self, *args):
        """ """
        assert np.shape(args)[0] == self.ndim, \
            " <CoordinateTransformation> : We need {} inputs in {}D.".format(self.ndim, self.ndim)
        assert all([np.ndim(args[i]) == 1 for i in range(self.ndim)]), \
            " <CoordinateTransformation> : We need 1-d inputs, will be meshgrided later."
        for i in range(self.ndim):
            assert np.max(args[i]) <= 1 and np.min(args[i]) >= -1, \
                " <CoordinateTransformation> : points_grid[{}]={} is wrong.".format(i, args[i])
        xietasigma = np.meshgrid(*args, indexing='ij')       
        self._set_evaluation_points_(*xietasigma)
        self._evaluation_points_grid_ = args
        
    def _set_evaluation_points_(self, *args):
        assert np.shape(args)[0] == self.ndim, \
            " <CoordinateTransformation> : evaluation points wrong."
        shape_0 = np.shape(args[0])
        assert all([np.shape(args[i]) == shape_0 for i in range(1, self.ndim)]), \
            " <CoordinateTransformation> : evaluation points wrong."
        for i in range(self.ndim):
            assert np.max(args[i]) <= 1 and np.min(args[i]) >= -1, \
                " <CoordinateTransformation> : points[{}] is wrong.".format(i)
        self._evaluation_points_ = args
    #----------------------------------------------------------------------------------
    
    @property
    def mapping(self):
        """ 
        Here we return the mapping. Since mapping maybe used at multiple places, we 
        save it to _mapping_.
        
        """
        assert self.evaluation_points is not None, \
            " <CoordinateTransformation> : no evaluation_points"
        return self.___compute_mapping___(None)
    
    def ___compute_mapping___(self, evaluation_points, element=None):
        """"""
        if evaluation_points is None:
            evaluation_points = self.evaluation_points
        if isinstance(element, int):
            return self.___compute_mapping_element_i___(evaluation_points, element)
        elif element is None:
            xyz = [np.zeros((self._mesh_.elements.global_num, *np.shape(evaluation_points[j])))
                        for j in range(self.ndim)]
            for i in range(self._mesh_.elements.global_num):
                xyz_i = self.___compute_mapping_element_i___(evaluation_points, i)
                for j in range(self.ndim): 
                    xyz[j][i] = xyz_i[j]
            return xyz
        else:
            raise Exception()
    
    def ___compute_mapping_element_i___(self, evaluation_points, i):
        """ """
        if evaluation_points is None:
            evaluation_points = self.evaluation_points
        region_name, local_indices = \
            self._mesh_.do.find.region_name_and_local_indices_of_element(i)
        origin, delta = \
            self._mesh_.do.find.reference_origin_and_size_of_element_of_given_local_indices(
                region_name, local_indices)
        xyz_i = self._mesh_.domain.regions(region_name).interpolation(
                *[(evaluation_points[j]+1)*0.5*delta[j] + origin[j] for j in range(self.ndim)])
        return xyz_i
    
    @property
    def Jacobian_matrix(self):
        """
        Here we return the Jacobian matrix. As Jacobian matrix may be used 
        for multiple times, we save it at _Jacobian_matrix_.
        
        """
        assert self.evaluation_points is not None, " <CoordinateTransformation> : no evaluation_points"
        return self.___compute_Jacobian_matrix___(self.evaluation_points)

    def ___compute_Jacobian_matrix___(self, evaluation_points, element=None):
        """
        Here we compute the Jacobian matrix.
        
        Parameters
        ----------
        evaluation_points : tuple
            A tuplde of shape = (self.ndim, ). At which points we evaluate the mapping
        element : NoneType or int
            If None (default), do it for all elements. If int, do it for that element.
            
        Returns
        -------
        XYZ_xietasigma : tuple
            A tuple of shape (ndim, ndim), and XYZ_xietasigma[0][0][k] represents 
            dx_dxi(*args) in kth element. and it is of shape (*shape(args[0])).
            
            So in 2D:
                _Jacobian_matrix_ = [[dx_dxi, dx_deta], 
                                     [dy_dxi, dy_deta]].
                
            While in 3D:
                _Jacobian_matrix_ = [[dx_dxi, dx_deta, dx_dsigma], 
                                     [dy_dxi, dy_deta, dy_dsigma], 
                                     [dz_dxi, dz_deta, dz_dsigma]].
        
        """
        if isinstance(element, int):
            return self.___compute_Jacobian_matrix_i___(evaluation_points, element)
        elif element is None:
            XYZ_xietasigma = [[np.zeros((self._mesh_.elements.global_num, *np.shape(evaluation_points[j])))
                                for _ in range(self.ndim)] for j in range(self.ndim)]
            for i in range(self._mesh_.elements.global_num):
                xyz_xietasigma = self.___compute_Jacobian_matrix_i___(evaluation_points, i)
                for j in range(self.ndim): 
                    for l in range(self.ndim):
                        XYZ_xietasigma[j][l][i] = xyz_xietasigma[j][l]
            return XYZ_xietasigma
        else:
            raise Exception
    
    def ___compute_Jacobian_matrix_i___(self, evaluation_points, element):
        """Compute the Jacobian matrix for ith element."""
        xyz_xietasigma = [[np.zeros(np.shape(evaluation_points[j])) 
                            for _ in range(self.ndim)]
                                for j in range(self.ndim)]
        i = element
        
        region_name, local_indices = \
            self._mesh_.do.find.region_name_and_local_indices_of_element(i)
        origin, delta = \
            self._mesh_.do.find.reference_origin_and_size_of_element_of_given_local_indices(
                region_name, local_indices)
        xyz_rst = self._mesh_.domain.regions(region_name).interpolation.Jacobian_matrix(
                *[(evaluation_points[j]+1)*0.5*delta[j] + origin[j] for j in range(self.ndim)])
        for j in range(self.ndim): 
            for l in range(self.ndim):
                xyz_xietasigma[j][l] = np.array([delta[l]/2])*xyz_rst[j][l]
        return xyz_xietasigma
        
    @property
    def metric(self):
        """
        The metric. Since our Jacobian and inverse of Jacobian are both square, we know
        that the metric g is equal to square of det(J). g = (det(J))**2 is due to the 
        fact that the Jacobian matrix is square. The definition of g usually is given 
        as g:= det(G) where G is the metric matrix, or metric tensor.
        
        """
        # noinspection PyUnresolvedReferences
        return self.Jacobian**2
    
    @property
    def metric_matrix(self):
        """
        Also called matric tensor. Let J be the Jacobian matrix. The metric_matrix is 
        denoted by G, G := J^T.dot(J). It is different from the metric:
            g := (det(J))**2 or g := det(G). 
        For a square Jacobian matrix, the metric turns out to be the square of the 
        determinant of the Jacobian matrix.
        
        The entries of G is normally denoted as g_{i,j}.
        
        g_{i,j}.
        
        """
        J = self.Jacobian_matrix
        G = [[None for _ in range(self.ndim)] for _ in range(self.ndim)]
        for i in range(self.ndim):
            for j in range(i, self.ndim):
                G[i][j] = J[0][i] * J[0][j]
                for l in range(1, self.ndim):
                    G[i][j] += J[l][i] * J[l][j]
                if i != j:
                    G[j][i] = G[i][j]
        return G

    @property
    def inverse_metric_matrix(self):
        """
        The inverse_metric_matrix is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.
        
        The entries of G^{-1} is normally denoted as g^{i,j}.
        
        g^{i,j}.
        
        """
        # noinspection PyUnresolvedReferences
        iJ = self.inverse_Jacobian_matrix
        iG = [[None for _ in range(self.ndim)] for _ in range(self.ndim)]
        for i in range(self.ndim):
            for j in range(i, self.ndim):
                iG[i][j] = iJ[i][0] * iJ[j][0]
                for l in range(1, self.ndim):
                    iG[i][j] += iJ[i][l] * iJ[j][l]
                if i != j:
                    iG[j][i] = iG[i][j]
        return iG