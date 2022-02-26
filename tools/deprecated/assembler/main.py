# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('/')
from importlib import import_module
from screws.frozen import FrozenOnly




class GatheringMatrix(FrozenOnly):
    def __init_subclass__(cls):
        pass

    @property
    def IS_gathering_matrix(self):
        return True






class Assembler(FrozenOnly):
    """
    Assembler.

    :param Gi_:
    :param G_j:
    """
    def __init__(self, Gi_, G_j=None):
        if G_j is None: G_j = Gi_
        self.___PRIVATE_check_Gi__G_j___(Gi_, G_j)
        self._assembler_ = None
        self._freeze_self_()

    def __call__(self, M, method=None, **kwargs):
        self._assembler_ = self.___PRIVATE_obtain_assembler___(M)
        if method is None:
            method =self._assembler_.default_method
        self._assembler_.___PRIVATE_check_M___(M)
        MA = getattr(self._assembler_, 'DO_assembling_van_'+method)(M, **kwargs)
        return MA

    def ___PRIVATE_check_Gi__G_j___(self, Gi_, G_j):
        if hasattr(Gi_, 'standard_properties') and '3dCSCG_form' in Gi_.standard_properties.tags:
            Gi_ = Gi_.numbering.gathering
        else:
            pass
        if hasattr(G_j, 'standard_properties') and '3dCSCG_form' in G_j.standard_properties.tags:
            G_j = G_j.numbering.gathering
        else:
            pass
        assert Gi_.IS_gathering_matrix, "Gi_ is not a gathering matrix object."
        assert G_j.IS_gathering_matrix, "G_j is not a gathering matrix object."
        self._gtypes_ = (Gi_.__class__.__name__, G_j.__class__.__name__)
        self._Gi__ = Gi_
        self._G_j_ = G_j

    def ___PRIVATE_obtain_assembler___(self, M):
        if M.__class__.__name__ == 'EWC_SparseMatrix':
            assembler_name = '_3dCSCG_EWC'
        else:
            raise Exception()

        if self._assembler_ is None or assembler_name != self._assembler_.__class__.__name__:
            return getattr(import_module('TOOLS.assembler.assemblers'), assembler_name)(self)
        else:
            return self._assembler_

    @property
    def Gi_(self):
        return self._Gi__

    @property
    def G_j(self):
        return self._G_j_




if __name__ == '__main__':
    # mpiexec python TOOLS\assembler\main.py
    import numpy as np
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller
    mesh = MeshGenerator('crazy', c=0.1)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z) + t/2

    scalar = FC('scalar', p)
    f0 = FC('0-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)

    f0.TW.func.body = scalar
    f0.TW.___DO_push_all_to_instant___()
    f0.discretize()

    M = f0.matrices.mass
    W = f3.operators.wedge(f0)

    M = Assembler(f0, f0)(M)
    W = Assembler(f0, f3)(W)
