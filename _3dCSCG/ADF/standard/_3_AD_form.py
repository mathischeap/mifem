




import sys
if './' not in sys.path: sys.path.append('./')


from _3dCSCG.ADF.standard.base.main import _3dCSCG_Algebra_DUAL_Standard_Form


class _3dCSCG_S3_ADF(_3dCSCG_Algebra_DUAL_Standard_Form):
    """
    Standard a dual 3-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        if name is None: name = orientation + '-oriented-3-ADF'
        super().__init__(3, mesh, space, orientation, name)
        self._prime_ = prime
        self._k_ = 3
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_3form')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()




    def ___PRIVATE_reset_cache___(self):
        """"""
        super().___PRIVATE_reset_cache___()




if __name__ == "__main__":

    # mpiexec -n 6 python _3dCSCG\ADF\standard\_3_AD_form.py
    import numpy as np
    from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy', c=0.0)([10, 10, 10], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([3,3,3], show_info=True)


    def p(t, x, y, z):
        return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def u(t, x, y, z):
        return 2*np.pi*np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def v(t, x, y, z):
        return 2*np.pi*np.sin(2*np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z) + t
    def w(t, x, y, z):
        return 2*np.pi*np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(2*np.pi*z) + t


    FC = FormCaller(mesh, space)
    df3 = FC('3-adf', numbering_parameters={'scheme_name': 'Naive', })
    dt2 = FC('2-adt', numbering_parameters={'scheme_name': 'Naive', })

    scalar = FC('scalar', p)
    vector = FC('vector', (u, v, w))

    df3.prime.TW.func.do.set_func_body_as(scalar)
    df3.prime.TW.current_time = 0
    df3.prime.TW.do.push_all_to_instant()
    df3.prime.do.discretize()

    dt2.prime.TW.func.do.set_func_body_as(scalar)
    dt2.prime.TW.current_time = 0
    dt2.prime.TW.do.push_all_to_instant()
    dt2.prime.do.discretize()


    df2 = df3.coboundary(dt2)

    df2.prime.TW.func.do.set_func_body_as(vector)
    df2.prime.TW.current_time = 0
    df2.prime.TW.do.push_all_to_instant()




    print(df3.prime.error.L())
    print(df2.prime.error.L())

