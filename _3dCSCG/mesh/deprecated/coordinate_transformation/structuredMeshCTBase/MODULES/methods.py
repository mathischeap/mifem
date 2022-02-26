# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Wed May 22 21:13:56 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
class CTMODMethods:
    """ """

    # noinspection PyUnresolvedReferences
    def ___inverse_metric_matrix_method___(self, iJ):
        """
        The inverse_matric_matrix is the matric matrix of the inverse Jacobian
        matrix or the metric of the inverse mapping. It is usally denoted as
        G^{-1}.
        
        The entries of G^{-1} is normally denoted as g^{i,j}.
        
        g^{i,j}.
        
        Parameters
        ----------
        iJ :
            The inverse_Jacobian_matrix.
        
        """
        iG = [[None for _ in range(self.ndim)] for _ in range(self.ndim)]
        for i in range(self.ndim):
            for j in range(i, self.ndim):
                iG[i][j] = iJ[i][0] * iJ[j][0]
                for l in range(1, self.ndim):
                    iG[i][j] += iJ[i][l] * iJ[j][l]
                if i != j:
                    iG[j][i] = iG[i][j]
        return iG
    