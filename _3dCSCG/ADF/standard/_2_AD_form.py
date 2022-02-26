"""

"""


import sys
if './' not in sys.path: sys.path.append('./')



from _3dCSCG.ADF.standard.base import _3dCSCG_Algebra_DUAL_Standard_Form






class _3dCSCG_S2_ADF(_3dCSCG_Algebra_DUAL_Standard_Form):
    """
    Standard a dual 2-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        if name is None: name = orientation + '-oriented-2-ADF'
        super().__init__(3, mesh, space, orientation, name)
        self._prime_ = prime
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_2form')
        self.RESET_cache()
        self._freeze_self_()




    def RESET_cache(self):
        """"""
        super().RESET_cache()









if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\_2_AD_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy',c=0.0, bounds=([-1,1],[-1,1],[-1,1]))(
                                        [12, 12, 12], EDM=None, show_info=True)

    space = SpaceInvoker('polynomials')([4, 4, 4], show_info=True)

    FC = FormCaller(mesh, space)