

"""
Sub-geometries of a mesh. Like we can pick a point, a slice or a volume from a mesh. Such a sub-geometry should
consist of sub-geometries of one or several elements.

"""

from root.config import *
from screws.frozen import FrozenOnly



class _3dCSCG_Mesh_SubGeometry(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh

        self._freeze_self_()



    def GENERATE_perpendicular_slice_object(self, *args, **kwargs):
        return _3dCSCG_MeshPerpendicularSlice(self._mesh_, *args, **kwargs)







class _3dCSCG_MeshPerpendicularSlice(FrozenOnly):
    """
    A particular slice-sub-geometry.

    It is a slice of a mesh. Such a sub-geometry should consist of slices of one or several elements.
    """
    def __init__(self, mesh, inputs, *args, **kwargs):
        """"""
        self._mesh_ = mesh

        if inputs.__class__.__name__ == 'RegionPerpendicularSlice':
            self._element_slices_, self._PTA_ = \
                self.___PRIVATE_generate_element_slices_from_A_RegionPerpendicularSlice__(inputs)
        else:
            raise NotImplementedError(f"Currently, we can only generate MeshSlice instance from RegionSlice. "
                                      f"No use {args} and {kwargs}")

        self._freeze_self_()


    def ___PRIVATE_generate_element_slices_from_A_RegionPerpendicularSlice__(self, RS):
        """

        :param RS: A RegionSlice instance.
        :return:
        """
        element_slices = dict()

        rn = RS._region_.name

        assert rn in self._mesh_.domain.regions, "regions does not exist."
        assert self._mesh_.domain.regions[rn] is RS._region_, "regions does not match."

        r, t, s = RS.r, RS.t, RS.s
        spacing = self._mesh_.elements.spacing[rn]

        A = self.___PRIVATE_locate_index_in_spacing___(r, spacing[0])
        B = self.___PRIVATE_locate_index_in_spacing___(t, spacing[1])
        C = self.___PRIVATE_locate_index_in_spacing___(s, spacing[2])

        _ = (A, B, C)
        PTA = 'rts'.index(RS.perpendicular_to_axis)
        assert len(_[PTA]) == 1, "something is wrong."

        PTA_position = (r, t, s)[PTA]

        IND = list()

        PTA_info = None

        NUM = list()
        for k, abc in enumerate(_):
            if len(abc) == 1:
                loc = abc[0][0]
                LS = spacing[k][loc:loc+2]
                start, end = LS
                assert start <= PTA_position < end, f"[start, end] = [{start}, {end}], PTA_position={PTA_position} wrong."
                PTA_element_position = -1 + 2 * (PTA_position - start) / (end - start)

                IND.append(range(0, len(spacing[k])))
                PTA_info = (loc, PTA_element_position)

            elif len(abc) == 2:
                A1, A2 = abc
                IND.append(range(A1[0], A2[1]))

                NUM.append(A2[1] - A1[0])

            else:
                raise Exception()

        ind1, ind2, ind3 = IND

        # Find local involved elements from IND ------------- BELOW -------------------------------------------

        involved_loc_elements = list()
        for e in self._mesh_.elements:
            r_n, ijk = self._mesh_.DO.FIND_region_name_and_local_indices_of_element(e)
            if r_n != rn:
                pass
            else:
                i, j, k = ijk
                if i in ind1 and j in ind2 and k in ind3:
                    pta = ijk[PTA]
                    if pta == PTA_info[0]:
                        involved_loc_elements.append(e)
                    else:
                        pass
                else:
                    pass

        num_loc = len(involved_loc_elements)
        ALL_num_loc = cOmm.gather(num_loc, root=mAster_rank)
        if rAnk == mAster_rank: assert sum(ALL_num_loc) == np.prod(NUM), "NUM checking failed."

        position = PTA_info[1]

        for e in involved_loc_elements:
            element = self._mesh_.elements[e]
            if PTA == 0:
                ESG_PS = element.sub_geometry.GENERATE_perpendicular_slice_object(xi = position)
            elif PTA == 1:
                ESG_PS = element.sub_geometry.GENERATE_perpendicular_slice_object(eta = position)
            elif PTA == 2:
                ESG_PS = element.sub_geometry.GENERATE_perpendicular_slice_object(sigma = position)
            else:
                raise Exception()

            element_slices[e] = ESG_PS

        return element_slices, ['xi', 'eta', 'sigma'][PTA]

    @staticmethod
    def ___PRIVATE_locate_index_in_spacing___(indices, spacing):
        """
        For example, given spacing = [0. ,   0.125, 0.25,  0.375, 0.5,   0.625, 0.75,  0.875, 1.],
        when indices = 0, then we return ([0,0],)

        when indices = 0.15, we return ([1,2],)

        when indices = [0.15, 0.8], we return ([1,2], [6,7],)

        :param indices:
        :param spacing:
        :return:
        """
        if isinstance(indices, (int, float)):
            indices = [indices,]

        location = tuple()
        for i in indices:
            LOC = None
            for j, s in enumerate(spacing):
                if i == s:
                    LOC = [j, j]
                    break
                elif s < i < spacing[j+1]:
                    LOC = [j, j+1]
                    break
                else:
                    pass
            assert LOC is not None, f"location={i} is wrong."
            location += (LOC,)

        return location


    def __getitem__(self, i):
        return self._element_slices_[i]

    def __contains__(self, i):
        return i in self._element_slices_

    def __iter__(self):
        for i in self._element_slices_:
            yield i

    def __len__(self):
        return len(self._element_slices_)

    @property
    def perpendicular_to_axis(self):
        return self._PTA_