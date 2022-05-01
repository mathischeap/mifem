
"""

"""
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.main import FrozenOnly
from screws.exceptions import LocalCochainShapeError
from root.config.main import *




class _3dCSCG_Edge_Cochain(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._local_ = None
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._local_EEW_ = None



    def ___PRIVATE_DO_gather_local_to_master___(self):
        """Do what the method name says."""
        local = cOmm.gather(self.local, root=mAster_rank)
        if rAnk == mAster_rank:
            LOCAL = dict()
            for li in local:
                if li is not None:
                    LOCAL.update(li)
            return LOCAL

    def ___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___(self):
        """make it regions-wise-element-local-indexed, thus we can save it and when read a form, we can always have the
        correct local cochain allocated even element numbering is different.
        """
        assert self.local is not None, "I have no local cochain!"

        RN_LI_dict = dict()
        for i in self.local:
            rn, loc_ind = self._ef_.mesh.do.find.region_name_and_local_indices_of_element(i)

            RN_LI_dict[i] = rn + '=|=' + str(loc_ind)

        RN_LI_dict = cOmm.gather(RN_LI_dict, root=mAster_rank)
        if rAnk == mAster_rank:
            RID = dict()
            for rid in RN_LI_dict:
                RID.update(rid)
        del RN_LI_dict

        LOCAL = self.___PRIVATE_DO_gather_local_to_master___()
        if rAnk == mAster_rank:

            RW_LOCAL = dict()

            for i in range(self._ef_.mesh.elements.GLOBAL_num):
                assert i in LOCAL, "something is wrong."
                # noinspection PyUnboundLocalVariable
                assert i in RID, "something is wrong."

                rn_loc_ind = RID[i]
                RW_LOCAL[rn_loc_ind] = LOCAL[i]

            return RW_LOCAL

    def ___PRIVATE_do_distribute_region_wise_local_index_grouped_cochain_to_local___(self, RW_LI_COCHAIN):
        """When we have the Region-wised local index grouped cochain, we can use this method to distribute it to local
        cochain. The regions-wise grouped cochain must be a full cochain in all cores.

        :param RW_LI_COCHAIN:
        :return:
        """

        LOC_COCHAIN = dict()

        for i in self._ef_.mesh.elements:
            rn, loc_ind = self._ef_.mesh.do.find.region_name_and_local_indices_of_element(i)
            dict_key = rn + '=|=' + str(loc_ind)
            LOC_COCHAIN[i] = RW_LI_COCHAIN[dict_key]

        self.local = LOC_COCHAIN






    def __getitem__(self, item):
        return self.local[item]

    def __contains__(self, item):
        return item in self.local

    def __iter__(self):
        for i in self.local:
            yield i

    def __len__(self):
        return len(self.local)

    #-- DEPENDENT PROPERTIES (MAJOR): When set local, clear BRANCHES by set branches to None ---------------------
    @property
    def local(self):
        """
        The local cochain.

        :return: None or A dict whose keys are local mesh element indices
            and values are cochain in corresponding mesh-elements.
        :rtype: Dict[int, numpy.ndarray]
        """
        return self._local_

    @local.setter
    def local(self, local):
        numOfElements = self._ef_.mesh.elements.num
        numOfBasis = self._ef_.num.basis
        try:
            assert isinstance(local, dict)
            assert len(local) == numOfElements
            for i in self._ef_.mesh.elements:
                assert np.shape(local[i]) == (numOfBasis,)
        except AssertionError:
            raise LocalCochainShapeError()

        self.___PRIVATE_reset_cache___()
        self._local_ = local

    #-- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set below, update local ------
    @property
    def local_EEW(self):
        """
        The local cochain in edge element.

        EEW stands for edge-Element-Wise.

        :return: A dict whose keys are edge element names and values are cochain in corresponding edge elements.
        :rtype: Dict[str, numpy.ndarray]
        """
        # this is important: do not use ``local_EEW`` or ``local``.
        if self._local_EEW_ is None and self._local_ is not None:
            self.___local_2_local_EEW___()
        return self._local_EEW_

    @local_EEW.setter
    def local_EEW(self, local_EEW):
        """
        :param local_EEW:
        :return:
        """
        NUM_basis_components = self._ef_.num.basis_components
        try:
            assert isinstance(local_EEW, dict)
            for key in local_EEW:
                ee = self._ef_.mesh.edge.elements[key]
                Cce = ee.CHARACTERISTIC_corner_edge
                assert local_EEW[key].shape == NUM_basis_components[Cce]
        except AssertionError:
            raise LocalCochainShapeError()

        self.___PRIVATE_reset_cache___()
        self._local_EEW_ = local_EEW
        self.___local_EEW_2_local___()


    def ___local_EEW_2_local___(self):
        """"""
        MAP = self._ef_.mesh.edge.elements.map
        for i in MAP:
            for key in MAP[i]:
                assert key in self._local_EEW_, "'local_TEW' is not full."

        local = dict()
        for i in MAP:
            local[i] = list()
            for key in MAP[i]:
                local[i].append(self._local_EEW_[key])
            local[i] = np.concatenate(local[i])
            assert local[i].shape == (self._ef_.num.basis,)
        self._local_ = local



    def ___local_2_local_EEW___(self):
        """"""
        raise NotImplementedError()






if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\cochain.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    scalar = FC('scalar', p)

    e0.TW.func.do.set_func_body_as(scalar)
    e0.TW.current_time = 0
    e0.TW.do.push_all_to_instant()

    e0.discretize()