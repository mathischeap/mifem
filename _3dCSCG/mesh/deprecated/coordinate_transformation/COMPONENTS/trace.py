# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Fri May  3 09:44:11 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from _3dCSCG.mesh.deprecated.coordinate_transformation.structuredMeshCTBase.COMPONENTS.trace import CTEXTTBase
from _3dCSCG.mesh.deprecated.coordinate_transformation.SEC import CoordinateTransformationSingleElementComputer


class CoordinateTransformationTrace(CTEXTTBase):
    """ 
    Here we store the mapping, Jacobian, inverse Jacobian, metric and inverse metric of 
    the trace (faces) of 3D mesh.
    
    """
    def __init__(self, ct):
        """ 
        Parameters
        ----------
        ct : CoordinateTransformation3D
        
        """
        super().__init__(ct)
    
    def ___generate_trace_evaluation_points___(self):
        """
        To do trace coordinate transformation, we need `evaluation_points_grid`
        rather than `evaluation_points`.
        
        Returns
        -------
        tep_N, tep_S, tep_W, tep_E, tep_B, tep_F : 
            trace evaluation points @North, South, West, East, Back and Front.
        
        """
        assert self._ct_.evaluation_points_grid is not None, \
            " <CoordinateTransformation3D.Trace3D> : no evaluation_points_grid"
        epg = self._ct_.evaluation_points_grid
        tep_N = [item.ravel('F') for item in np.meshgrid([-1], epg[1], epg[2], indexing='ij')]
        tep_S = [item.ravel('F') for item in np.meshgrid([+1], epg[1], epg[2], indexing='ij')]
        tep_W = [item.ravel('F') for item in np.meshgrid(epg[0], [-1], epg[2], indexing='ij')]
        tep_E = [item.ravel('F') for item in np.meshgrid(epg[0], [+1], epg[2], indexing='ij')]
        tep_B = [item.ravel('F') for item in np.meshgrid(epg[0], epg[1], [-1], indexing='ij')]
        tep_F = [item.ravel('F') for item in np.meshgrid(epg[0], epg[1], [+1], indexing='ij')]
        return tep_N, tep_S, tep_W, tep_E, tep_B, tep_F
    
    @property
    def mapping(self):
        """ 
        The mapping for 3D trace elements (2D geometry objects in 3D.). 
        
        Returns
        -------
        self._mapping_ : dict
            keys: trace element numbering; values: mapping coordinates. The 
            cooridnates are in ndarray of shape (3, ...) corresponding to x, y, 
            z.
            
        """
#        if self._mapping_ is None:
#            tep_N, tep_S, tep_W, tep_E, tep_B, tep_F = self.___generate_trace_evaluation_points___()
#            M = {'N':self._ct_.___compute_mapping___(tep_N), # North, 
#                 'S':self._ct_.___compute_mapping___(tep_S), # South, 
#                 'W':self._ct_.___compute_mapping___(tep_W), # West, 
#                 'E':self._ct_.___compute_mapping___(tep_E), # East, 
#                 'B':self._ct_.___compute_mapping___(tep_B), # Back, 
#                 'F':self._ct_.___compute_mapping___(tep_F)} # Front
#            self._mapping_ = {}
#            for i in self.mesh.trace.elements.position_representive:
#                # Notice that here we actually will only use the info in 
#                # element_position_representive, that means a lot of data in 
#                # `M` will not be used. But I believe this is not a big problem
#                # as the computing data amount of the trace mapping is much 
#                # less  than computing data of the mesh elements mapping. So we
#                # leave it as it is now.
#                element, side = self.mesh.trace.elements.position_representive[i].split('-')
#                element = int(element)
#                m = M[side]
#                self._mapping_[i] = np.array([m[0][element, ...], m[1][element, ...], m[2][element, ...]])
#        return self._mapping_
        tep_N, tep_S, tep_W, tep_E, tep_B, tep_F = self.___generate_trace_evaluation_points___()
        M = {'N':self._ct_.___compute_mapping___(tep_N), # North, 
             'S':self._ct_.___compute_mapping___(tep_S), # South, 
             'W':self._ct_.___compute_mapping___(tep_W), # West, 
             'E':self._ct_.___compute_mapping___(tep_E), # East, 
             'B':self._ct_.___compute_mapping___(tep_B), # Back, 
             'F':self._ct_.___compute_mapping___(tep_F)} # Front
        _mapping_ = {}
        for i in self._mesh_.trace.elements.position_representive:
            # Notice that here we actually will only use the info in 
            # element_position_representive, that means a lot of data in 
            # `M` will not be used. But I believe this is not a big problem
            # as the computing data amount of the trace mapping is much 
            # less  than computing data of the mesh elements mapping. So we
            # leave it as it is now.
            element, side = self._mesh_.trace.elements.position_representive[i].split('-')
            element = int(element)
            m = M[side]
            _mapping_[i] = np.array([m[0][element, ...], m[1][element, ...], m[2][element, ...]])
            
        return _mapping_
    
    @property
    def Jacobian_matrix(self):
        """ 
        The Jacobian matrix for 3D trace elements (2D geometry objects in 3D 
        space).
        
        Returns
        -------
        _Jacobian_matrix_ : dict
            Keys: # trace element. Values: Jacobian_matrix of shape (3,2).
        
        """
#        if self._Jacobian_matrix_ is None:
#            tep_N, tep_S, tep_W, tep_E, tep_B, tep_F = self.___generate_trace_evaluation_points___()
#            J = {'N':self._ct_.___compute_Jacobian_matrix___(tep_N), # North, 
#                 'S':self._ct_.___compute_Jacobian_matrix___(tep_S), # South, 
#                 'W':self._ct_.___compute_Jacobian_matrix___(tep_W), # West, 
#                 'E':self._ct_.___compute_Jacobian_matrix___(tep_E), # East, 
#                 'B':self._ct_.___compute_Jacobian_matrix___(tep_B), # Back, 
#                 'F':self._ct_.___compute_Jacobian_matrix___(tep_F)} # Front
#            self._Jacobian_matrix_ = {}
#            for i in self.mesh.trace.elements.position_representive:
#                element, side = self.mesh.trace.elements.position_representive[i].split('-')
#                element = int(element)
#                j = J[side]
#                if side in 'NS':
#                    self._Jacobian_matrix_[i] = np.array(((j[0][1][element, ...], j[0][2][element, ...]),
#                                                          (j[1][1][element, ...], j[1][2][element, ...]),
#                                                          (j[2][1][element, ...], j[2][2][element, ...])))
#                elif side in 'WE':
#                    self._Jacobian_matrix_[i] = np.array(((j[0][0][element, ...], j[0][2][element, ...]),
#                                                          (j[1][0][element, ...], j[1][2][element, ...]),
#                                                          (j[2][0][element, ...], j[2][2][element, ...])))
#                elif side in 'BF':
#                    self._Jacobian_matrix_[i] = np.array(((j[0][0][element, ...], j[0][1][element, ...]),
#                                                          (j[1][0][element, ...], j[1][1][element, ...]),
#                                                          (j[2][0][element, ...], j[2][1][element, ...])))
#                else:
#                    raise Exception
#        return self._Jacobian_matrix_
        tep_N, tep_S, tep_W, tep_E, tep_B, tep_F = self.___generate_trace_evaluation_points___()
        J = {'N':self._ct_.___compute_Jacobian_matrix___(tep_N), # North, 
             'S':self._ct_.___compute_Jacobian_matrix___(tep_S), # South, 
             'W':self._ct_.___compute_Jacobian_matrix___(tep_W), # West, 
             'E':self._ct_.___compute_Jacobian_matrix___(tep_E), # East, 
             'B':self._ct_.___compute_Jacobian_matrix___(tep_B), # Back, 
             'F':self._ct_.___compute_Jacobian_matrix___(tep_F)} # Front
        _Jacobian_matrix_ = {}
        for i in self._mesh_.trace.elements.position_representive:
            element, side = self._mesh_.trace.elements.position_representive[i].split('-')
            element = int(element)
            j = J[side]
            if side in 'NS':
                _Jacobian_matrix_[i] = np.array(((j[0][1][element, ...], j[0][2][element, ...]),
                                                 (j[1][1][element, ...], j[1][2][element, ...]),
                                                 (j[2][1][element, ...], j[2][2][element, ...])))
            elif side in 'WE':
                _Jacobian_matrix_[i] = np.array(((j[0][0][element, ...], j[0][2][element, ...]),
                                                 (j[1][0][element, ...], j[1][2][element, ...]),
                                                 (j[2][0][element, ...], j[2][2][element, ...])))
            elif side in 'BF':
                _Jacobian_matrix_[i] = np.array(((j[0][0][element, ...], j[0][1][element, ...]),
                                                 (j[1][0][element, ...], j[1][1][element, ...]),
                                                 (j[2][0][element, ...], j[2][1][element, ...])))
            else:
                raise Exception
        return _Jacobian_matrix_
        
    @property
    def metric(self):
        """ g, which should be det(metric_matrix): det(G). """
        G = self.metric_matrix
        _metric_ = {}
        for i in self._mesh_.trace.elements.position_representive:
            _metric_[i] = G[i][0][0]*G[i][1][1] - G[i][0][1]*G[i][1][0]
        return _metric_    
        
    @property
    def inverse_metric_matrix(self):
        """ 
        In the trace case, it will be better to compute it by inverting the
        'metric_matrix'.
        
        And this only happens in 3D. Because in 2D, the trace is a 1-d stuff, 
        we will only needs to use metric g.
        
        """
        G = self.metric_matrix
        _imm_ = {}
        for i in self._mesh_.trace.elements.position_representive:
            _imm_[i] = np.linalg.inv(np.array(G[i]).transpose((2,0,1))).transpose((1,2,0))
        return _imm_
    
    @property
    def inverse_Jacobian_matrix(self):
        """ 
        The inverse Jacobian matrix for 3D trace elements (2D geometry objects 
        in 3D space).
        
        Returns
        -------
        iJm : dict
            Keys are the numbering of trace-elements.
        
        """
        tep_N, tep_S, tep_W, tep_E, tep_B, tep_F = self.___generate_trace_evaluation_points___()
        tep_ = {'N': tep_N, # North, 
                'S': tep_S, # South, 
                'W': tep_W, # West, 
                'E': tep_E, # East, 
                'B': tep_B, # Back, 
                'F': tep_F} # Front,
        SEC = CoordinateTransformationSingleElementComputer(self._ct_)
        smtepr = self._mesh_.trace.elements.position_representive
        tiJm = {}
        strNS, strWE, strBF = 'NS', 'WE', 'BF'
        for i in range(self._mesh_.trace.elements.num): 
            # `i`th trace-element 
            j, side = smtepr[i].split('-')
            j = int(j)
            iJ = SEC.___inverse_Jacobian_matrix_ep___(j, tep_[side])
            if side in strNS:
                tiJm[i] = np.array(((iJ[1][0], iJ[1][1], iJ[1][2]),
                                    (iJ[2][0], iJ[2][1], iJ[2][2])))
            elif side in strWE:
                tiJm[i] = np.array(((iJ[0][0], iJ[0][1], iJ[0][2]),
                                    (iJ[2][0], iJ[2][1], iJ[2][2])))
            elif side in strBF:
                tiJm[i] = np.array(((iJ[0][0], iJ[0][1], iJ[0][2]),
                                    (iJ[1][0], iJ[1][1], iJ[1][2])))
            else:
                raise Exception()
        return tiJm