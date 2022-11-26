# -*- coding: utf-8 -*-
"""
INTRO

<unittest> <unittests_Dump> <test_No2_1DumpForm>.

Yi Zhang (C)
Created on Fri May  3 09:42:55 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from components.freeze.main import FrozenOnly

# %% CoordinateTransformationDump
class CoordinateTransformationDump(FrozenOnly):
    """ 
    Here we have the coordinate transformations for the dump-elements.
    
    <unittest> <unittests_Dump> <test_No2_1DumpForm>.
    
    """
    def __init__(self, ct):
        """ 
        Parameters
        ----------
        ct : CoordinateTransformation3D
        
        """
        self._ct_ = ct
        self._freeze_self_()

    @property
    def _mesh_(self):
        return self._ct_._mesh_
    
    @property
    def ndim(self):
        """ dump element is a n-2 dimensonal object. """
        return self._ct_.ndim - 2
        
    def ___generate_dump_evaluation_points___(self):
        """
        To do dump coordinate transformation, we need `evaluation_points_grid` rather 
        than `evaluation_points`.
        
        <unittest> <unittests_Dump> <test_No2_1DumpForm>.
        
        Returns
        -------
        dict :
            A dict of keys: WB, EB, WF, EF, NB, SB, NF, SF, NW, SW, NE, SE: dump 
            evaluation points @West-Back, East-Back, ......
        
        """
        assert self._ct_.evaluation_points_grid is not None, \
            " <CoordinateTransformation3D.Trace3D> : no evaluation_points_grid"
        epg = self._ct_.evaluation_points_grid
        # ind = {0:'WB', 1:'EB',  2:'WF',  3:'EF',  
        #        4:'NB', 5:'SB',  6:'NF',  7:'SF', 
        #        8:'NW', 9:'SW', 10:'NE', 11:'SE'}
        WB = [item.ravel('F') for item in np.meshgrid(epg[0], [-1], [-1], indexing='ij')]
        EB = [item.ravel('F') for item in np.meshgrid(epg[0], [ 1], [-1], indexing='ij')]
        WF = [item.ravel('F') for item in np.meshgrid(epg[0], [-1], [ 1], indexing='ij')]
        EF = [item.ravel('F') for item in np.meshgrid(epg[0], [ 1], [ 1], indexing='ij')]
        NB = [item.ravel('F') for item in np.meshgrid([-1], epg[1], [-1], indexing='ij')]
        SB = [item.ravel('F') for item in np.meshgrid([ 1], epg[1], [-1], indexing='ij')]
        NF = [item.ravel('F') for item in np.meshgrid([-1], epg[1], [ 1], indexing='ij')]
        SF = [item.ravel('F') for item in np.meshgrid([ 1], epg[1], [ 1], indexing='ij')]
        NW = [item.ravel('F') for item in np.meshgrid([-1], [-1], epg[2], indexing='ij')]
        SW = [item.ravel('F') for item in np.meshgrid([ 1], [-1], epg[2], indexing='ij')]
        NE = [item.ravel('F') for item in np.meshgrid([-1], [ 1], epg[2], indexing='ij')]
        SE = [item.ravel('F') for item in np.meshgrid([ 1], [ 1], epg[2], indexing='ij')]
        return {'WB':WB, 'EB':EB, 'WF':WF, 'EF':EF, 
                'NB':NB, 'SB':SB, 'NF':NF, 'SF':SF, 
                'NW':NW, 'SW':SW, 'NE':NE, 'SE':SE}
    
    @property
    def mapping(self):
        """
        <unittest> <unittests_Dump> <test_No2_1DumpForm>.
        
        """
        dep = self.___generate_dump_evaluation_points___()
        M = {'WB': self._ct_.___compute_mapping___(dep['WB']),
             'EB': self._ct_.___compute_mapping___(dep['EB']), 
             'WF': self._ct_.___compute_mapping___(dep['WF']), 
             'EF': self._ct_.___compute_mapping___(dep['EF']),  
             'NB': self._ct_.___compute_mapping___(dep['NB']), 
             'SB': self._ct_.___compute_mapping___(dep['SB']), 
             'NF': self._ct_.___compute_mapping___(dep['NF']), 
             'SF': self._ct_.___compute_mapping___(dep['SF']), 
             'NW': self._ct_.___compute_mapping___(dep['NW']), 
             'SW': self._ct_.___compute_mapping___(dep['SW']), 
             'NE': self._ct_.___compute_mapping___(dep['NE']),  
             'SE': self._ct_.___compute_mapping___(dep['SE'])}
        _mapping_ = {}
        depr = self._mesh_.dump.elements.position_representive
        for i in depr:
            # Notice that here we actually will only use the info in 
            # position_representive, that means a lot of data in `M` will not be used. 
            # But I believe this is not a big problem as the computing data amount of 
            # the dump mapping is much less than computing data of the mesh elements 
            # mapping. So we leave it as it is now.
            element, edge = depr[i].split('-')
            element = int(element)
            m = M[edge]
            _mapping_[i] = np.array([m[0][element, ...], m[1][element, ...], m[2][element, ...]])
        return _mapping_
    
    @property
    def Jacobian_matrix(self):
        """ 
        <unittest> <unittests_Dump> <test_No2_1DumpForm>.
        
        Returns
        -------
        _Jacobian_matrix_ : dict
            Keys: # trace element. Values: Jacobian_matrix of shape (3,1).
                
        """
        dep = self.___generate_dump_evaluation_points___()
        J = {'WB': self._ct_.___compute_Jacobian_matrix___(dep['WB']),
             'EB': self._ct_.___compute_Jacobian_matrix___(dep['EB']), 
             'WF': self._ct_.___compute_Jacobian_matrix___(dep['WF']), 
             'EF': self._ct_.___compute_Jacobian_matrix___(dep['EF']),  
             'NB': self._ct_.___compute_Jacobian_matrix___(dep['NB']), 
             'SB': self._ct_.___compute_Jacobian_matrix___(dep['SB']), 
             'NF': self._ct_.___compute_Jacobian_matrix___(dep['NF']), 
             'SF': self._ct_.___compute_Jacobian_matrix___(dep['SF']), 
             'NW': self._ct_.___compute_Jacobian_matrix___(dep['NW']), 
             'SW': self._ct_.___compute_Jacobian_matrix___(dep['SW']), 
             'NE': self._ct_.___compute_Jacobian_matrix___(dep['NE']),  
             'SE': self._ct_.___compute_Jacobian_matrix___(dep['SE'])}
        _Jacobian_matrix_ = {}
        depr = self._mesh_.dump.elements.position_representive
        for i in depr:
            # Notice that here we actually will only use the info in 
            # position_representive, that means a lot of data in `M` will not be used. 
            # But I believe this is not a big problem as the computing data amount of 
            # the dump mapping is much less than computing data of the mesh elements 
            # mapping. So we leave it as it is now.
            element, edge = depr[i].split('-')
            element = int(element)
            j = J[edge]
            if edge in ('WB', 'EB', 'WF', 'EF'): # dx
                _Jacobian_matrix_[i] = np.array(((j[0][0][element, ...],),
                                                 (j[1][0][element, ...],),
                                                 (j[2][0][element, ...],)))
            elif edge in ('NB', 'SB', 'NF', 'SF'): # dy
                _Jacobian_matrix_[i] = np.array(((j[0][1][element, ...],),
                                                 (j[1][1][element, ...],),
                                                 (j[2][1][element, ...],)))
            elif edge in ('NW', 'SW', 'NE', 'SE'): # dz
                _Jacobian_matrix_[i] = np.array(((j[0][2][element, ...],),
                                                 (j[1][2][element, ...],),
                                                 (j[2][2][element, ...],)))
            else:
                raise Exception()
        return _Jacobian_matrix_
    
    @property
    def metric_matrix(self):
        """ 
        The entries of metric_matrix is normally denoted as g_{i,j}.  For the 
        dump-element, it should be of shape (1, 1).
        
        G[k][0][0] refers to the the metrix_matrix. Since the matric_matrix is of shape
        (1,1). So, to refer the metric_matrix, we do G[k][0][0], k is `k`th 
        dump-element.
        
        <unittest> <unittests_Dump> <test_No2_1DumpForm>.
        
        """
        J = self.Jacobian_matrix
        G = {}
        for k in self._mesh_.dump.elements.position_representive:
            G[k] = J[k][0][0]*J[k][0][0] + J[k][1][0]*J[k][1][0] + J[k][2][0]*J[k][2][0]
            G[k] = ((G[k],),)
        return G
    
    @property
    def metric(self):
        """ 
        the metric g.
        
        <unittest> <unittests_Dump> <test_No2_1DumpForm>.
        
        """
        G = self.metric_matrix
        g = {}
        for k in self._mesh_.dump.elements.position_representive:
            g[k] = G[k][0][0]
        return g
    
    @property
    def inverse_metric_matrix(self):
        """ The inverse metric matrix."""
        G = self.metric_matrix
        iG = {}
        for k in self._mesh_.dump.elements.position_representive:
            iG[k] = np.reciprocal(G[k])
        return iG