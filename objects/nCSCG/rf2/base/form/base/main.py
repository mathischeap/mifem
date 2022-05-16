# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 6:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass
from objects.nCSCG.rf2.base.form.base.TW.main import nCSCG_RF2_FormTW


class nCSCG_RF2_FormBase(FrozenClass):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._k_ = None
        self.___is_wrapped_in_ADF___ = False # change this when we wrap it into an ADF
        self.___define_parameters___ = None

        self._TW_ = None
        self._func_ = None

    @property
    def mesh(self):
        return self._mesh_

    @property
    def k(self):
        return self._k_

    @property
    def ndim(self):
        """Return the dimensions."""
        return self.mesh.ndim

    #---------- must have properties -------------------------------------------
    @property
    def TW(self):
        if self._TW_ is None:
             self._TW_ = nCSCG_RF2_FormTW(self)
        return self._TW_

    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self):
        raise NotImplementedError()

    @property
    def cochain(self):
        raise NotImplementedError()

    #---------------------------------------------------------------------------------------
    @property
    def ___personal_parameters___(self):
        """
        Personal parameters are parameters that are additional to define_parameters. The full parameters
        are combination of both personal and define parameters.

        IMPORTANT:
        """
        if self.TW.func is None:
            func_parameters = None
        else:
            func_parameters = self.TW.func.___parameters___

        RGW_cochain_local = self.cochain.___Pr_RGW_cochain___()
        # RGW_cochain_local will only be in master, in slaves, it is None.
        return {'TW_current_time': self.TW.current_time, # to save a object, it must have current time.
                'TW_func_parameters': func_parameters,
                'region_wise_cochain_local': RGW_cochain_local}

    @property
    def ___parameters___(self):
        """
        ___parameters___ will be called when you call standard_properties.parameters.

        :return:
        """
        parameters = dict()
        parameters.update(self.___define_parameters___)
        parameters.update(self.___personal_parameters___)
        return parameters




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
