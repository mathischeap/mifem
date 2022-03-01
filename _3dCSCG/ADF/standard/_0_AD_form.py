
import sys
if './' not in sys.path: sys.path.append('./')


from _3dCSCG.ADF.standard.base import _3dCSCG_Algebra_DUAL_Standard_Form


class _3dCSCG_S0_ADF(_3dCSCG_Algebra_DUAL_Standard_Form):
    """
    Standard a dual 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        if name is None: name = orientation + '-oriented-0-ADF'
        super().__init__(3, mesh, space, orientation, name)
        self._prime_ = prime
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_0form')
        self.RESET_cache()
        self._freeze_self_()



    def RESET_cache(self):
        """"""
        super().RESET_cache()









if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\standard\_0_AD_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',2)], show_info=True)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    FC = FormCaller(mesh, space)
    df0 = FC('0-adf', numbering_parameters={'scheme_name': 'Naive', }, name='potential')

    df0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    df0.prime.do.discretize()

