# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Wed May 22 21:13:56 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
class CTMODGenerators:
    """ """
    # noinspection PyUnresolvedReferences
    def ___mapping_generator___(self, evaluation_points):
        """ """
        k = 0
        while k < self._mesh_.elements.num:
            yield self.___compute_mapping_element_i___(evaluation_points, k)
            k += 1

    # noinspection PyUnresolvedReferences
    def ___Jacobian_matrix_generator___(self, evaluation_points):
        """ """
        k = 0
        while k < self._mesh_.elements.num:
            yield self.___compute_Jacobian_matrix_i___(evaluation_points, k)
            k += 1

    # noinspection PyUnresolvedReferences
    def ___inverse_metric_matrix_generator___(self, evaluation_points):
        """ """
        iJg = self.___inverse_Jacobian_matrix_generator___(evaluation_points)
        k = 0
        while k < self._mesh_.elements.num:
            iJ = iJg.__next__()
            iG = [[None for _ in range(self.ndim)] for _ in range(self.ndim)]
            for i in range(self.ndim):
                for j in range(i, self.ndim):
                    iG[i][j] = iJ[i][0] * iJ[j][0]
                    for l in range(1, self.ndim):
                        iG[i][j] += iJ[i][l] * iJ[j][l]
                    if i != j:
                        iG[j][i] = iG[i][j]
            yield iG
            k += 1
        