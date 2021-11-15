# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Sat May  4 14:07:51 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from SCREWS.frozen import FrozenOnly
from scipy.sparse import csc_matrix

class IncidenceMatrix(FrozenOnly):
    """ 
    Clearly, the incidence matrix of a form only depends on the local numbering and the basis function degree. We have
    fixed the way of numbering local dofs. That is the reason why we can evaluate basis here. Therefore, we can already
    compute the incidence matrix.
    
    """
    def __init__(self, FS):
        """ """
        assert FS.ndim == 3, " <IncidenceMatrix> "
        self._FS_ = FS
        self._freeze_self_()
    
    @property
    def _0Form(self):
        """
        Here we generate the incidence matrix for 0-form in 3D.
        """
        sn = self._FS_.local_numbering._0Form
        dn = self._FS_.local_numbering._1Form
        E = np.zeros((self._FS_.num_basis._1Form[0], self._FS_.num_basis._0Form[0]), dtype=int)
        
        I, J, K = np.shape(dn[0])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[0][i,j,k], sn[0][i,j,k]]   = -1   # North
                    E[dn[0][i,j,k], sn[0][i+1,j,k]] = +1   # South
        
        I, J, K = np.shape(dn[1])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[1][i,j,k], sn[0][i,j,k]]   = -1    # West
                    E[dn[1][i,j,k], sn[0][i,j+1,k]] = +1    # East
        
        I, J, K = np.shape(dn[2])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[2][i,j,k], sn[0][i,j,k]]   = -1    # Back
                    E[dn[2][i,j,k], sn[0][i,j,k+1]] = +1    # Front
        
        return csc_matrix(E)
    
    @property
    def _1Form(self):
        """
        Here we generate the incidence matrix for 1-form in 3D.
        """
        sn = self._FS_.local_numbering._1Form
        dn = self._FS_.local_numbering._2Form
        E = np.zeros((self._FS_.num_basis._2Form[0], self._FS_.num_basis._1Form[0]), dtype=int)
        
        I, J, K = np.shape(dn[0])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[0][i,j,k], sn[1][i,j,k  ]] = +1   # Back
                    E[dn[0][i,j,k], sn[1][i,j,k+1]] = -1   # Front
                    E[dn[0][i,j,k], sn[2][i,j  ,k]] = -1   # West
                    E[dn[0][i,j,k], sn[2][i,j+1,k]] = +1   # East
        
        I, J, K = np.shape(dn[1])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[1][i,j,k], sn[0][i,j,k  ]] = -1    # Back
                    E[dn[1][i,j,k], sn[0][i,j,k+1]] = +1    # Front
                    E[dn[1][i,j,k], sn[2][i  ,j,k]] = +1    # North
                    E[dn[1][i,j,k], sn[2][i+1,j,k]] = -1    # South
        
        I, J, K = np.shape(dn[2])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[2][i,j,k], sn[0][i,j  ,k]] = +1    # West
                    E[dn[2][i,j,k], sn[0][i,j+1,k]] = -1    # East
                    E[dn[2][i,j,k], sn[1][i  ,j,k]] = -1    # North
                    E[dn[2][i,j,k], sn[1][i+1,j,k]] = +1    # South
        return csc_matrix(E)
    
    @property
    def _2Form(self):
        """
        Here we generate the incidence matrix for 2-form in 2D.
        """
        sn = self._FS_.local_numbering._2Form
        dn = self._FS_.local_numbering._3Form
        E = np.zeros((self._FS_.num_basis._3Form[0], self._FS_.num_basis._2Form[0]), dtype=int)
        
        I, J, K = np.shape(dn[0])
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    E[dn[0][i,j,k], sn[0][i  ,j,k]] = -1 # North
                    E[dn[0][i,j,k], sn[0][i+1,j,k]] = +1 # South
                    E[dn[0][i,j,k], sn[1][i,j  ,k]] = -1 # West
                    E[dn[0][i,j,k], sn[1][i,j+1,k]] = +1 # East
                    E[dn[0][i,j,k], sn[2][i,j,k  ]] = -1 # Back
                    E[dn[0][i,j,k], sn[2][i,j,k+1]] = +1 # Front
        return csc_matrix(E)