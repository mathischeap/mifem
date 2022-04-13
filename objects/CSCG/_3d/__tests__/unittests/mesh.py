# -*- coding: utf-8 -*-
"""
Mesh related unittests.
"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from screws.quadrature import Quadrature
from screws.exceptions import ThreeDimensionalTransfiniteInterpolationError
from objects.CSCG._3d.mesh.domain.inputs.allocator import DomainInputAllocator
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
import random
import os
from objects.CSCG._3d.__tests__.random_objects.form_caller import random_mesh_of_elements_around


def test_Mesh_NO0_element_division_and_numbering_quality():
    """"""
    if rAnk == mAster_rank:
        print("~~~ [test_Mesh_NO0_element_division_and_numbering_quality] ...... ", flush=True)

    try:
        MESH = MeshGenerator('LDC', l=1, w=1.2, h=1.5)([f'Lobatto:{13}', f'Lobatto:{14}', f'Lobatto:{15}'], EDM='debug')
        mesh = MeshGenerator('LDC', l=1, w=1.2, h=1.5)([f'Lobatto:{13}', f'Lobatto:{14}', f'Lobatto:{15}'])

        if 6 <= sIze <= 24:
            A = MESH.___PRIVATE_element_division_and_numbering_quality___()[0]
            B = mesh.___PRIVATE_element_division_and_numbering_quality___()[0]
            assert A <= B, "Smarter division should result in better quality."

    except ThreeDimensionalTransfiniteInterpolationError:

        if rAnk == mAster_rank:
            print("   = Partial test SKIPPED.", flush=True)

    MESH = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='debug')
    if sIze >= 4:
        mesh = MeshGenerator('bridge_arch_cracked',)([3,2,4], EDM='SWV0')
    else:
        mesh = MeshGenerator('bridge_arch_cracked', )([3, 2, 4])
    A = MESH.___PRIVATE_element_division_and_numbering_quality___()[0]
    B = mesh.___PRIVATE_element_division_and_numbering_quality___()[0]

    if sIze <= 24:
        assert A <= B, "Smarter division should result in better quality."


    return 1


def test_Mesh_NO1_mesh_general():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO1_mesh_general} ...... ", flush=True)

    # test method ___PRIVATE_do_find_slave_of_element___ ...
    mesh = MeshGenerator('crazy')([5, 4, 3], EDM='debug')
    for i in range(mesh.elements.GLOBAL_num):
        sn = mesh.do.find.slave_of_element(i)
        assert i in mesh._element_distribution_[sn]
    mesh = MeshGenerator('crazy')([1, 2, 1], EDM='debug')
    for i in range(mesh.elements.GLOBAL_num):
        sn = mesh.do.find.slave_of_element(i)
        assert i in mesh._element_distribution_[sn]
    return 1


def test_Mesh_NO2_trace_elements():
    """Unittests for the trace elements."""
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO2_trace_elements} ...... ", flush=True)

    mesh = MeshGenerator('crazy')([2, 2, 2], EDM='debug')
    trace_elements = mesh.trace.elements
    benchmark = {0: ('0N', 'North'), 1: ('0S', '1N'), 2: ('0W', 'West'), 3: ('0E', '2W'),
                 4: ('0B', 'Back'), 5: ('0F', '4B'), 6: ('1S', 'South'), 7: ('1W', 'West'),
                 8: ('1E', '3W'), 9: ('1B', 'Back'), 10: ('1F', '5B'), 11: ('2N', 'North'),
                 12: ('2S', '3N'), 13: ('2E', 'East'), 14: ('2B', 'Back'), 15: ('2F', '6B'),
                 16: ('3S', 'South'), 17: ('3E', 'East'), 18: ('3B', 'Back'), 19: ('3F', '7B'),
                 20: ('4N', 'North'), 21: ('4S', '5N'), 22: ('4W', 'West'), 23: ('4E', '6W'),
                 24: ('4F', 'Front'), 25: ('5S', 'South'), 26: ('5W', 'West'), 27: ('5E', '7W'),
                 28: ('5F', 'Front'), 29: ('6N', 'North'), 30: ('6S', '7N'), 31: ('6E', 'East'),
                 32: ('6F', 'Front'), 33: ('7S', 'South'), 34: ('7E', 'East'), 35: ('7F', 'Front')}
    for i in trace_elements:
        tei = trace_elements[i]
        assert tei.positions == benchmark[i], \
            f"trace element [{i}] position {tei.positions} != benchmark {benchmark[i]}"

    benchmark = {6: [29, 30, 23, 31, 15, 32], 7: [30, 33, 27, 34, 19, 35],
                 2: [11, 12, 3, 13, 14, 15], 3: [12, 16, 8, 17, 18, 19],
                 4: [20, 21, 22, 23, 5, 24], 5: [21, 25, 26, 27, 10, 28],
                 0: [0, 1, 2, 3, 4, 5], 1: [1, 6, 7, 8, 9, 10]}

    for i in trace_elements.map:
        assert trace_elements.map[i] == benchmark[i]

    benchmark = {0: ('0N', '1S'), 1: ('0S', '1N'), 2: ('0W', '2E'), 3: ('0E', '2W'), 4: ('0B', '4F'),
                 5: ('0F', '4B'), 6: ('1W', '3E'), 7: ('1E', '3W'), 8: ('1B', '5F'), 9: ('1F', '5B'),
                 10: ('2N', '3S'), 11: ('2S', '3N'), 12: ('2B', '6F'), 13: ('2F', '6B'),
                 14: ('3B', '7F'), 15: ('3F', '7B'), 16: ('4N', '5S'), 17: ('4S', '5N'),
                 18: ('4W', '6E'), 19: ('4E', '6W'), 20: ('5W', '7E'), 21: ('5E', '7W'),
                 22: ('6N', '7S'), 23: ('6S', '7N')}

    mesh = MeshGenerator('crazy_periodic')([2, 2, 2], EDM='debug')
    trace_elements = mesh.trace.elements
    for i in trace_elements:
        tei = trace_elements[i]
        assert tei.positions == benchmark[i]

    benchmark = {0: [0, 1, 2, 3, 4, 5], 1: [1, 0, 6, 7, 8, 9], 2: [10, 11, 3, 2, 12, 13],
                 3: [11, 10, 7, 6, 14, 15], 4: [16, 17, 18, 19, 5, 4], 5: [17, 16, 20, 21, 9, 8],
                 6: [22, 23, 19, 18, 13, 12], 7: [23, 22, 21, 20, 15, 14]}
    for i in trace_elements.map:
        assert trace_elements.map[i] == benchmark[i]
    benchmark = {0 : ('0N', '1S') ,
                 1 : ('0S', '1N') ,
                 2 : ('0W', '0E') ,
                 3 : ('0B', '0F') ,
                 4 : ('1W', '1E') ,
                 5 : ('1B', '1F')}
    mesh = MeshGenerator('crazy_periodic')([2, 1, 1], EDM='debug')
    trace_elements = mesh.trace.elements
    for i in trace_elements:
        tei = trace_elements[i]
        assert tei.positions == benchmark[i]
    benchmark = {0: [0, 1, 2, 2, 3, 3], 1: [1, 0, 4, 4, 5, 5]}
    for i in trace_elements.map:
        assert trace_elements.map[i] == benchmark[i]

    return 1


def test_Mesh_NO2a_trace_elements_CT():
    """Unittests for the coordinate transformation of trace elements."""
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO2a_trace_elements_CT} ...... ", flush=True)

    if rAnk == mAster_rank:
        while 1:
            el1 = random.randint(2,5)
            el2 = random.randint(2,5)
            el3 = random.randint(2,5)
            if el1 * el2 * el3 < 100: # do not test too big mesh
                break
        c = random.uniform(0.0, 0.3)
        if c < 0.15:
            c = 0
        _i_ = random.randint(3,6)
        _j_ = random.randint(3,6)
    else:
        el1, el2, el3, c, _i_, _j_ = [None for _ in range(6)]
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)
    _i_, _j_ = cOmm.bcast([_i_, _j_], root=mAster_rank)
    m = MeshGenerator('crazy', c=c)([el1, el2, el3], EDM='debug')

    xi = np.linspace(-1, 1, _i_)
    et = np.linspace(-1, 1, _j_)
    xi, et = np.meshgrid(xi, et, indexing='ij')
    tes = m.trace.elements
    POSITION = dict()
    MAPPING = dict()
    METRIC = dict()
    JM = dict()
    MM = dict()
    for i in tes:
        te = tes[i]
        if te.IS.shared_by_cores:
            POSITION[i] = te.CHARACTERISTIC_side
            MAPPING[i] = te.coordinate_transformation.mapping(xi, et)
            METRIC[i] = te.coordinate_transformation.metric(xi, et)
            JM[i] = te.coordinate_transformation.Jacobian_matrix(xi, et)
            iJM = te.coordinate_transformation.inverse_Jacobian_matrix(xi, et)
            MM[i] = te.coordinate_transformation.metric_matrix(xi, et)

            J00, J01 = JM[i][0]
            J10, J11 = JM[i][1]
            J20, J21 = JM[i][2]

            if J00.__class__.__name__ == 'ndarray':
                assert J00.shape == (_i_, _j_)
            if J10.__class__.__name__ == 'ndarray':
                assert J10.shape == (_i_, _j_)
            if J20.__class__.__name__ == 'ndarray':
                assert J20.shape == (_i_, _j_)
            if J01.__class__.__name__ == 'ndarray':
                assert J01.shape == (_i_, _j_)
            if J11.__class__.__name__ == 'ndarray':
                assert J11.shape == (_i_, _j_)
            if J21.__class__.__name__ == 'ndarray':
                assert J21.shape == (_i_, _j_)

            iJ00, iJ01, iJ02 = iJM[0]
            iJ10, iJ11, iJ12 = iJM[1]
            iJJ00 = iJ00*J00 + iJ01*J10 + iJ02*J20
            np.testing.assert_array_almost_equal(iJJ00, 1)
            iJJ11 = iJ10*J01 + iJ11*J11 + iJ12*J21
            np.testing.assert_array_almost_equal(iJJ11, 1)
            iJJ01 = iJ00*J01 + iJ01*J11 + iJ02*J21
            iJJ10 = iJ10*J00 + iJ11*J10 + iJ12*J20
            np.testing.assert_array_almost_equal(iJJ01, 0)
            np.testing.assert_array_almost_equal(iJJ10, 0)

        else:
            jm = te.coordinate_transformation.Jacobian_matrix(xi, et)
            ijm = te.coordinate_transformation.inverse_Jacobian_matrix(xi, et)

            J00, J01 = jm[0]
            J10, J11 = jm[1]
            J20, J21 = jm[2]

            if J00.__class__.__name__ == 'ndarray':
                assert J00.shape == (_i_, _j_)
            if J10.__class__.__name__ == 'ndarray':
                assert J10.shape == (_i_, _j_)
            if J20.__class__.__name__ == 'ndarray':
                assert J20.shape == (_i_, _j_)
            if J01.__class__.__name__ == 'ndarray':
                assert J01.shape == (_i_, _j_)
            if J11.__class__.__name__ == 'ndarray':
                assert J11.shape == (_i_, _j_)
            if J21.__class__.__name__ == 'ndarray':
                assert J21.shape == (_i_, _j_)

            iJ00, iJ01, iJ02 = ijm[0]
            iJ10, iJ11, iJ12 = ijm[1]
            iJJ00 = iJ00*J00 + iJ01*J10 + iJ02*J20
            np.testing.assert_array_almost_equal(iJJ00, 1)
            iJJ11 = iJ10*J01 + iJ11*J11 + iJ12*J21
            np.testing.assert_array_almost_equal(iJJ11, 1)
            iJJ01 = iJ00*J01 + iJ01*J11 + iJ02*J21
            iJJ10 = iJ10*J00 + iJ11*J10 + iJ12*J20
            np.testing.assert_array_almost_equal(iJJ01, 0)
            np.testing.assert_array_almost_equal(iJJ10, 0)

    POSITION = cOmm.gather(POSITION, root=mAster_rank)
    MAPPING = cOmm.gather(MAPPING, root=sEcretary_rank)
    METRIC = cOmm.gather(METRIC, root=sEcretary_rank)
    MM = cOmm.gather(MM, root=mAster_rank)
    JM = cOmm.gather(JM, root=mAster_rank)

    if rAnk == mAster_rank: #to check we get same results in different cores.
        _POS_ = dict()
        for PI in POSITION:
            for i in PI:
                if i in _POS_:
                    _POS_[i] += PI[i]
                else:
                    _POS_[i] = PI[i]
        check_tuple = ('NS', 'SN', 'WE', 'EW', 'FB', 'BF')
        for i in _POS_:
            assert _POS_[i] in check_tuple, \
                f"trace element No. [{i}] position wrong."


        _MM_ = dict()
        for MI in MM:
            for i in MI:
                if i in _MM_:
                    _MM_[i] += (MI[i],)
                else:
                    _MM_[i] = (MI[i],)

        _JM_ = dict()
        for MI in JM:
            for i in MI:
                if i in _JM_:
                    _JM_[i] += (MI[i],)
                else:
                    _JM_[i] = (MI[i],)

        for i in _MM_:
            assert len(_MM_[i]) == 2
            assert len(_JM_[i]) == 2
            # noinspection PyTupleAssignmentBalance
            A, B = _MM_[i]
            _00_01_, _10_11_ = A
            a00, a01 = _00_01_
            a10, a11 = _10_11_
            _00_01_, _10_11_ = B
            b00, b01 = _00_01_
            b10, b11 = _10_11_
            np.testing.assert_almost_equal( np.sum(np.abs(a00-b00)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a01-b01)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a10-b10)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a11-b11)), 0)

            # noinspection PyTupleAssignmentBalance
            A, B = _JM_[i]
            _0_, _1_, _2_ = A
            a00, a01 = _0_
            a10, a11 = _1_
            a20, a21 = _2_
            _0_, _1_, _2_ = B
            b00, b01 = _0_
            b10, b11 = _1_
            b20, b21 = _2_
            np.testing.assert_almost_equal( np.sum(np.abs(a00-b00)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a01-b01)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a10-b10)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a11-b11)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a20-b20)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(a21-b21)), 0)

    if rAnk == sEcretary_rank: #to check we get same results in different cores.
        _MAP_ = dict()
        for MI in MAPPING:
            for i in MI:
                if i in _MAP_:
                    _MAP_[i] += (MI[i],)
                else:
                    _MAP_[i] = (MI[i],)

        for i in _MAP_:
            assert len(_MAP_[i]) == 2
            # noinspection PyTupleAssignmentBalance
            A, B = _MAP_[i]
            x, y, z = A
            a, b, c = B
            np.testing.assert_almost_equal( np.sum(np.abs(x-a)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(y-b)), 0)
            np.testing.assert_almost_equal( np.sum(np.abs(z-c)), 0)

        _MET_ = dict()
        for MI in METRIC:
            for i in MI:
                if i in _MET_:
                    _MET_[i] += (MI[i],)
                else:
                    _MET_[i] = (MI[i],)
        for i in _MET_:
            assert len(_MET_[i]) == 2
            # noinspection PyTupleAssignmentBalance
            A, B = _MET_[i]
            np.testing.assert_almost_equal(np.sum(np.abs(A - B)), 0)

    xi = np.linspace(-1, 1, _i_ + 1)
    et = np.linspace(-1, 1, _j_ + 2)
    xi, et = np.meshgrid(xi, et, indexing='ij')
    for i in tes:
        te = tes[i]
        jm = te.coordinate_transformation.Jacobian_matrix(xi, et)
        ijm = te.coordinate_transformation.inverse_Jacobian_matrix(xi, et)

        J00, J01 = jm[0]
        J10, J11 = jm[1]
        J20, J21 = jm[2]

        if J00.__class__.__name__ == 'ndarray':
            assert J00.shape == (_i_+1, _j_+2)
        if J10.__class__.__name__ == 'ndarray':
            assert J10.shape == (_i_+1, _j_+2)
        if J20.__class__.__name__ == 'ndarray':
            assert J20.shape == (_i_+1, _j_+2)
        if J01.__class__.__name__ == 'ndarray':
            assert J01.shape == (_i_+1, _j_+2)
        if J11.__class__.__name__ == 'ndarray':
            assert J11.shape == (_i_+1, _j_+2)
        if J21.__class__.__name__ == 'ndarray':
            assert J21.shape == (_i_+1, _j_+2)

        iJ00, iJ01, iJ02 = ijm[0]
        iJ10, iJ11, iJ12 = ijm[1]
        iJJ00 = iJ00*J00 + iJ01*J10 + iJ02*J20
        np.testing.assert_array_almost_equal(iJJ00, 1)
        iJJ11 = iJ10*J01 + iJ11*J11 + iJ12*J21
        np.testing.assert_array_almost_equal(iJJ11, 1)
        iJJ01 = iJ00*J01 + iJ01*J11 + iJ02*J21
        iJJ10 = iJ10*J00 + iJ11*J10 + iJ12*J20
        np.testing.assert_array_almost_equal(iJJ01, 0)
        np.testing.assert_array_almost_equal(iJJ10, 0)

    return 1


def test_Mesh_NO3_elements_CT():
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO3_elements_CT} ...... ", flush=True)

    if rAnk == mAster_rank:
        el1 = random.randint(1,4)
        el2 = random.randint(1,3)
        el3 = random.randint(2,3)
        c = random.uniform(0, 0.3)
        if c < 0.15:
            c = 0
    else:
        el1, el2, el3, c = None, None, None, None
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)
    m = MeshGenerator('crazy_periodic', c=c)([el1, el2, el3], EDM='debug')
    m.___PRIVATE_generate_element_global_numbering___()

    if rAnk == mAster_rank:
        r = np.linspace(random.uniform(-1, -0.9), random.uniform(0.95, 0.99), random.randint(2,4))
        s = np.linspace(random.uniform(-1, -0.8), random.uniform(0.85, 0.9), random.randint(1,3))
        t = np.linspace(random.uniform(-1, -0.85), random.uniform(0.88, 0.93), random.randint(1,5))
    else:
        r, s, t = None, None, None

    r, s, t = cOmm.bcast([r, s, t], root=mAster_rank)

    r,s,t = np.meshgrid(r,s,t, indexing='ij')

    m.___TEST_MODE___ = True
    m.___DEPRECATED_ct___.evaluated_at(r, s, t)
    mapping = m.___DEPRECATED_ct___.mapping
    JM = m.___DEPRECATED_ct___.Jacobian_matrix
    J = m.___DEPRECATED_ct___.Jacobian
    iJM = m.___DEPRECATED_ct___.inverse_Jacobian_matrix
    iJ = m.___DEPRECATED_ct___.inverse_Jacobian

    M = m.___DEPRECATED_ct___.metric
    MM = m.___DEPRECATED_ct___.metric_matrix
    iMM = m.___DEPRECATED_ct___.inverse_metric_matrix


    _mapping = m.elements.coordinate_transformation.mapping(r, s, t)
    _X = m.elements.coordinate_transformation.X(r, s, t)
    _Y = m.elements.coordinate_transformation.Y(r, s, t)
    _Z = m.elements.coordinate_transformation.Z(r, s, t)
    _JM = m.elements.coordinate_transformation.Jacobian_matrix(r, s, t)
    _J00 = m.elements.coordinate_transformation.J00(r, s, t)
    _J01 = m.elements.coordinate_transformation.J01(r, s, t)
    _J02 = m.elements.coordinate_transformation.J02(r, s, t)
    _J10 = m.elements.coordinate_transformation.J10(r, s, t)
    _J11 = m.elements.coordinate_transformation.J11(r, s, t)
    _J12 = m.elements.coordinate_transformation.J12(r, s, t)
    _J20 = m.elements.coordinate_transformation.J20(r, s, t)
    _J21 = m.elements.coordinate_transformation.J21(r, s, t)
    _J22 = m.elements.coordinate_transformation.J22(r, s, t)
    _J = m.elements.coordinate_transformation.Jacobian(r, s, t, J=_JM)
    _M = m.elements.coordinate_transformation.metric(r, s, t, detJ=_J)
    _MM = m.elements.coordinate_transformation.metric_matrix(r, s, t, J=_JM)
    _iJM = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t, J=_JM)
    _iJ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t, iJ=_iJM)
    _iMM = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t, iJ=_iJM)

    for i in m.elements:
        ei = m.elements[i]

        mapping_i = ei.coordinate_transformation.mapping(r,s,t)

        X = ei.coordinate_transformation.X(r, s, t)
        Y = ei.coordinate_transformation.Y(r, s, t)
        Z = ei.coordinate_transformation.Z(r, s, t)

        np.testing.assert_array_almost_equal(mapping[0][i], X)
        np.testing.assert_array_almost_equal(mapping[1][i], Y)
        np.testing.assert_array_almost_equal(mapping[2][i], Z)

        np.testing.assert_array_almost_equal(mapping[0][i], mapping_i[0])
        np.testing.assert_array_almost_equal(mapping[1][i], mapping_i[1])
        np.testing.assert_array_almost_equal(mapping[2][i], mapping_i[2])

        JM_i = ei.coordinate_transformation.Jacobian_matrix(r,s,t)

        np.testing.assert_array_almost_equal(JM[0][0][i], JM_i[0][0])
        np.testing.assert_array_almost_equal(JM[0][1][i], JM_i[0][1])
        np.testing.assert_array_almost_equal(JM[0][2][i], JM_i[0][2])
        np.testing.assert_array_almost_equal(JM[1][0][i], JM_i[1][0])
        np.testing.assert_array_almost_equal(JM[1][1][i], JM_i[1][1])
        np.testing.assert_array_almost_equal(JM[1][2][i], JM_i[1][2])
        np.testing.assert_array_almost_equal(JM[2][0][i], JM_i[2][0])
        np.testing.assert_array_almost_equal(JM[2][1][i], JM_i[2][1])
        np.testing.assert_array_almost_equal(JM[2][2][i], JM_i[2][2])

        J00 = ei.coordinate_transformation.J00(r,s,t)
        J01 = ei.coordinate_transformation.J01(r,s,t)
        J02 = ei.coordinate_transformation.J02(r,s,t)
        J10 = ei.coordinate_transformation.J10(r,s,t)
        J11 = ei.coordinate_transformation.J11(r,s,t)
        J12 = ei.coordinate_transformation.J12(r,s,t)
        J20 = ei.coordinate_transformation.J20(r,s,t)
        J21 = ei.coordinate_transformation.J21(r,s,t)
        J22 = ei.coordinate_transformation.J22(r,s,t)

        np.testing.assert_array_almost_equal(JM[0][0][i], J00)
        np.testing.assert_array_almost_equal(JM[0][1][i], J01)
        np.testing.assert_array_almost_equal(JM[0][2][i], J02)
        np.testing.assert_array_almost_equal(JM[1][0][i], J10)
        np.testing.assert_array_almost_equal(JM[1][1][i], J11)
        np.testing.assert_array_almost_equal(JM[1][2][i], J12)
        np.testing.assert_array_almost_equal(JM[2][0][i], J20)
        np.testing.assert_array_almost_equal(JM[2][1][i], J21)
        np.testing.assert_array_almost_equal(JM[2][2][i], J22)

        J0 = ei.coordinate_transformation.J0_(r,s,t)
        J1 = ei.coordinate_transformation.J1_(r,s,t)
        J2 = ei.coordinate_transformation.J2_(r,s,t)
        np.testing.assert_array_almost_equal(J0[0], J00)
        np.testing.assert_array_almost_equal(J0[1], J01)
        np.testing.assert_array_almost_equal(J0[2], J02)
        np.testing.assert_array_almost_equal(J1[0], J10)
        np.testing.assert_array_almost_equal(J1[1], J11)
        np.testing.assert_array_almost_equal(J1[2], J12)
        np.testing.assert_array_almost_equal(J2[0], J20)
        np.testing.assert_array_almost_equal(J2[1], J21)
        np.testing.assert_array_almost_equal(J2[2], J22)

        J_i = ei.coordinate_transformation.Jacobian(r,s,t)
        iJ_i = ei.coordinate_transformation.inverse_Jacobian(r,s,t)
        M_i = ei.coordinate_transformation.metric(r,s,t)




        np.testing.assert_array_almost_equal(J[i], J_i)
        np.testing.assert_array_almost_equal(iJ[i], iJ_i)
        np.testing.assert_array_almost_equal(M[i], M_i)


        # test iJ @ J = I _________________________________________________
        iJM_i = ei.coordinate_transformation.inverse_Jacobian_matrix(r,s,t)
        iJ0, iJ1, iJ2 = iJM_i
        iJ00, iJ01, iJ02 = iJ0
        iJ10, iJ11, iJ12 = iJ1
        iJ20, iJ21, iJ22 = iJ2

        iJJ00 = iJ00*J00 + iJ01*J10 + iJ02*J20
        iJJ01 = iJ00*J01 + iJ01*J11 + iJ02*J21
        iJJ02 = iJ00*J02 + iJ01*J12 + iJ02*J22

        iJJ10 = iJ10*J00 + iJ11*J10 + iJ12*J20
        iJJ11 = iJ10*J01 + iJ11*J11 + iJ12*J21
        iJJ12 = iJ10*J02 + iJ11*J12 + iJ12*J22

        iJJ20 = iJ20*J00 + iJ21*J10 + iJ22*J20
        iJJ21 = iJ20*J01 + iJ21*J11 + iJ22*J21
        iJJ22 = iJ20*J02 + iJ21*J12 + iJ22*J22

        np.testing.assert_array_almost_equal(iJJ00, 1)
        np.testing.assert_array_almost_equal(iJJ01, 0)
        np.testing.assert_array_almost_equal(iJJ02, 0)

        np.testing.assert_array_almost_equal(iJJ10, 0)
        np.testing.assert_array_almost_equal(iJJ11, 1)
        np.testing.assert_array_almost_equal(iJJ12, 0)

        np.testing.assert_array_almost_equal(iJJ20, 0)
        np.testing.assert_array_almost_equal(iJJ21, 0)
        np.testing.assert_array_almost_equal(iJJ22, 1)
        #---------------------------------------------------------------

        np.testing.assert_array_almost_equal(iJM[0][0][i], iJM_i[0][0])
        np.testing.assert_array_almost_equal(iJM[0][1][i], iJM_i[0][1])
        np.testing.assert_array_almost_equal(iJM[0][2][i], iJM_i[0][2])
        np.testing.assert_array_almost_equal(iJM[1][0][i], iJM_i[1][0])
        np.testing.assert_array_almost_equal(iJM[1][1][i], iJM_i[1][1])
        np.testing.assert_array_almost_equal(iJM[1][2][i], iJM_i[1][2])
        np.testing.assert_array_almost_equal(iJM[2][0][i], iJM_i[2][0])
        np.testing.assert_array_almost_equal(iJM[2][1][i], iJM_i[2][1])
        np.testing.assert_array_almost_equal(iJM[2][2][i], iJM_i[2][2])

        MM_i = ei.coordinate_transformation.metric_matrix(r,s,t)
        iMM_i = ei.coordinate_transformation.inverse_metric_matrix(r,s,t)
        np.testing.assert_array_almost_equal(MM[0][0][i], MM_i[0][0])
        np.testing.assert_array_almost_equal(MM[0][1][i], MM_i[0][1])
        np.testing.assert_array_almost_equal(MM[0][2][i], MM_i[0][2])
        np.testing.assert_array_almost_equal(MM[1][0][i], MM_i[1][0])
        np.testing.assert_array_almost_equal(MM[1][1][i], MM_i[1][1])
        np.testing.assert_array_almost_equal(MM[1][2][i], MM_i[1][2])
        np.testing.assert_array_almost_equal(MM[2][0][i], MM_i[2][0])
        np.testing.assert_array_almost_equal(MM[2][1][i], MM_i[2][1])
        np.testing.assert_array_almost_equal(MM[2][2][i], MM_i[2][2])

        np.testing.assert_array_almost_equal(iMM[0][0][i], iMM_i[0][0])
        np.testing.assert_array_almost_equal(iMM[0][1][i], iMM_i[0][1])
        np.testing.assert_array_almost_equal(iMM[0][2][i], iMM_i[0][2])
        np.testing.assert_array_almost_equal(iMM[1][0][i], iMM_i[1][0])
        np.testing.assert_array_almost_equal(iMM[1][1][i], iMM_i[1][1])
        np.testing.assert_array_almost_equal(iMM[1][2][i], iMM_i[1][2])
        np.testing.assert_array_almost_equal(iMM[2][0][i], iMM_i[2][0])
        np.testing.assert_array_almost_equal(iMM[2][1][i], iMM_i[2][1])
        np.testing.assert_array_almost_equal(iMM[2][2][i], iMM_i[2][2])

        np.testing.assert_array_almost_equal(_mapping[i][0], mapping_i[0])
        np.testing.assert_array_almost_equal(_mapping[i][1], mapping_i[1])
        np.testing.assert_array_almost_equal(_mapping[i][2], mapping_i[2])
        np.testing.assert_array_almost_equal(_X[i], mapping_i[0])
        np.testing.assert_array_almost_equal(_Y[i], mapping_i[1])
        np.testing.assert_array_almost_equal(_Z[i], mapping_i[2])


        np.testing.assert_array_almost_equal(_JM[i][0][0], J00)
        np.testing.assert_array_almost_equal(_JM[i][0][1], J01)
        np.testing.assert_array_almost_equal(_JM[i][0][2], J02)
        np.testing.assert_array_almost_equal(_JM[i][1][0], J10)
        np.testing.assert_array_almost_equal(_JM[i][1][1], J11)
        np.testing.assert_array_almost_equal(_JM[i][1][2], J12)
        np.testing.assert_array_almost_equal(_JM[i][2][0], J20)
        np.testing.assert_array_almost_equal(_JM[i][2][1], J21)
        np.testing.assert_array_almost_equal(_JM[i][2][2], J22)

        np.testing.assert_array_almost_equal(_J00[i], J00)
        np.testing.assert_array_almost_equal(_J01[i], J01)
        np.testing.assert_array_almost_equal(_J02[i], J02)
        np.testing.assert_array_almost_equal(_J10[i], J10)
        np.testing.assert_array_almost_equal(_J11[i], J11)
        np.testing.assert_array_almost_equal(_J12[i], J12)
        np.testing.assert_array_almost_equal(_J20[i], J20)
        np.testing.assert_array_almost_equal(_J21[i], J21)
        np.testing.assert_array_almost_equal(_J22[i], J22)

        np.testing.assert_array_almost_equal(_J[i], J_i)
        np.testing.assert_array_almost_equal(_M[i], M_i)
        np.testing.assert_array_almost_equal(_iJ[i], iJ_i)

        np.testing.assert_array_almost_equal(_iJM[i][0][0], iJM_i[0][0])
        np.testing.assert_array_almost_equal(_iJM[i][0][1], iJM_i[0][1])
        np.testing.assert_array_almost_equal(_iJM[i][0][2], iJM_i[0][2])
        np.testing.assert_array_almost_equal(_iJM[i][1][0], iJM_i[1][0])
        np.testing.assert_array_almost_equal(_iJM[i][1][1], iJM_i[1][1])
        np.testing.assert_array_almost_equal(_iJM[i][1][2], iJM_i[1][2])
        np.testing.assert_array_almost_equal(_iJM[i][2][0], iJM_i[2][0])
        np.testing.assert_array_almost_equal(_iJM[i][2][1], iJM_i[2][1])
        np.testing.assert_array_almost_equal(_iJM[i][2][2], iJM_i[2][2])

        np.testing.assert_array_almost_equal(_MM[i][0][0], MM_i[0][0])
        np.testing.assert_array_almost_equal(_MM[i][0][1], MM_i[0][1])
        np.testing.assert_array_almost_equal(_MM[i][0][2], MM_i[0][2])
        np.testing.assert_array_almost_equal(_MM[i][1][0], MM_i[1][0])
        np.testing.assert_array_almost_equal(_MM[i][1][1], MM_i[1][1])
        np.testing.assert_array_almost_equal(_MM[i][1][2], MM_i[1][2])
        np.testing.assert_array_almost_equal(_MM[i][2][0], MM_i[2][0])
        np.testing.assert_array_almost_equal(_MM[i][2][1], MM_i[2][1])
        np.testing.assert_array_almost_equal(_MM[i][2][2], MM_i[2][2])

        np.testing.assert_array_almost_equal(_iMM[i][0][0], iMM_i[0][0])
        np.testing.assert_array_almost_equal(_iMM[i][0][1], iMM_i[0][1])
        np.testing.assert_array_almost_equal(_iMM[i][0][2], iMM_i[0][2])
        np.testing.assert_array_almost_equal(_iMM[i][1][0], iMM_i[1][0])
        np.testing.assert_array_almost_equal(_iMM[i][1][1], iMM_i[1][1])
        np.testing.assert_array_almost_equal(_iMM[i][1][2], iMM_i[1][2])
        np.testing.assert_array_almost_equal(_iMM[i][2][0], iMM_i[2][0])
        np.testing.assert_array_almost_equal(_iMM[i][2][1], iMM_i[2][1])
        np.testing.assert_array_almost_equal(_iMM[i][2][2], iMM_i[2][2])

    m = MeshGenerator('crazy_periodic', c=0.)(element_layout=[el1, el2, el3], EDM='debug')
    m.___TEST_MODE___ = True
    m.___PRIVATE_generate_element_global_numbering___()
    _mapping = m.elements.coordinate_transformation.mapping(r, s, t)
    _X = m.elements.coordinate_transformation.X(r, s, t)
    _Y = m.elements.coordinate_transformation.Y(r, s, t)
    _Z = m.elements.coordinate_transformation.Z(r, s, t)
    _JM = m.elements.coordinate_transformation.Jacobian_matrix(r, s, t)
    _J00 = m.elements.coordinate_transformation.J00(r, s, t)
    _J01 = m.elements.coordinate_transformation.J01(r, s, t)
    _J02 = m.elements.coordinate_transformation.J02(r, s, t)
    _J10 = m.elements.coordinate_transformation.J10(r, s, t)
    _J11 = m.elements.coordinate_transformation.J11(r, s, t)
    _J12 = m.elements.coordinate_transformation.J12(r, s, t)
    _J20 = m.elements.coordinate_transformation.J20(r, s, t)
    _J21 = m.elements.coordinate_transformation.J21(r, s, t)
    _J22 = m.elements.coordinate_transformation.J22(r, s, t)
    _J = m.elements.coordinate_transformation.Jacobian(r, s, t, J=_JM)
    _M = m.elements.coordinate_transformation.metric(r, s, t, detJ=_J)
    _MM = m.elements.coordinate_transformation.metric_matrix(r, s, t, J=_JM)
    _iJM = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t, J=_JM)
    _iJ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t, iJ=_iJM)
    _iMM = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t, iJ=_iJM)

    for i in m.elements:
        ei = m.elements[i]

        mapping_i = ei.coordinate_transformation.mapping(r,s,t)
        np.testing.assert_array_almost_equal(_mapping[i][0], mapping_i[0])
        np.testing.assert_array_almost_equal(_mapping[i][1], mapping_i[1])
        np.testing.assert_array_almost_equal(_mapping[i][2], mapping_i[2])
        np.testing.assert_array_almost_equal(_X[i], mapping_i[0])
        np.testing.assert_array_almost_equal(_Y[i], mapping_i[1])
        np.testing.assert_array_almost_equal(_Z[i], mapping_i[2])


        J00 = ei.coordinate_transformation.J00(r,s,t)
        J01 = ei.coordinate_transformation.J01(r,s,t)
        J02 = ei.coordinate_transformation.J02(r,s,t)
        J10 = ei.coordinate_transformation.J10(r,s,t)
        J11 = ei.coordinate_transformation.J11(r,s,t)
        J12 = ei.coordinate_transformation.J12(r,s,t)
        J20 = ei.coordinate_transformation.J20(r,s,t)
        J21 = ei.coordinate_transformation.J21(r,s,t)
        J22 = ei.coordinate_transformation.J22(r,s,t)
        np.testing.assert_array_almost_equal(_JM[i][0][0], J00)
        np.testing.assert_array_almost_equal(_JM[i][0][1], J01)
        np.testing.assert_array_almost_equal(_JM[i][0][2], J02)
        np.testing.assert_array_almost_equal(_JM[i][1][0], J10)
        np.testing.assert_array_almost_equal(_JM[i][1][1], J11)
        np.testing.assert_array_almost_equal(_JM[i][1][2], J12)
        np.testing.assert_array_almost_equal(_JM[i][2][0], J20)
        np.testing.assert_array_almost_equal(_JM[i][2][1], J21)
        np.testing.assert_array_almost_equal(_JM[i][2][2], J22)

        np.testing.assert_array_almost_equal(_J00[i], J00)
        np.testing.assert_array_almost_equal(_J01[i], J01)
        np.testing.assert_array_almost_equal(_J02[i], J02)
        np.testing.assert_array_almost_equal(_J10[i], J10)
        np.testing.assert_array_almost_equal(_J11[i], J11)
        np.testing.assert_array_almost_equal(_J12[i], J12)
        np.testing.assert_array_almost_equal(_J20[i], J20)
        np.testing.assert_array_almost_equal(_J21[i], J21)
        np.testing.assert_array_almost_equal(_J22[i], J22)

        J_i = ei.coordinate_transformation.Jacobian(r,s,t)
        iJ_i = ei.coordinate_transformation.inverse_Jacobian(r,s,t)
        M_i = ei.coordinate_transformation.metric(r,s,t)
        np.testing.assert_array_almost_equal(_J[i], J_i)
        np.testing.assert_array_almost_equal(_M[i], M_i)
        np.testing.assert_array_almost_equal(_iJ[i], iJ_i)

        iJM_i = ei.coordinate_transformation.inverse_Jacobian_matrix(r,s,t)
        MM_i = ei.coordinate_transformation.metric_matrix(r,s,t)
        iMM_i = ei.coordinate_transformation.inverse_metric_matrix(r,s,t)

        np.testing.assert_array_almost_equal(_iJM[i][0][0], iJM_i[0][0])
        np.testing.assert_array_almost_equal(_iJM[i][0][1], iJM_i[0][1])
        np.testing.assert_array_almost_equal(_iJM[i][0][2], iJM_i[0][2])
        np.testing.assert_array_almost_equal(_iJM[i][1][0], iJM_i[1][0])
        np.testing.assert_array_almost_equal(_iJM[i][1][1], iJM_i[1][1])
        np.testing.assert_array_almost_equal(_iJM[i][1][2], iJM_i[1][2])
        np.testing.assert_array_almost_equal(_iJM[i][2][0], iJM_i[2][0])
        np.testing.assert_array_almost_equal(_iJM[i][2][1], iJM_i[2][1])
        np.testing.assert_array_almost_equal(_iJM[i][2][2], iJM_i[2][2])

        np.testing.assert_array_almost_equal(_MM[i][0][0], MM_i[0][0])
        np.testing.assert_array_almost_equal(_MM[i][0][1], MM_i[0][1])
        np.testing.assert_array_almost_equal(_MM[i][0][2], MM_i[0][2])
        np.testing.assert_array_almost_equal(_MM[i][1][0], MM_i[1][0])
        np.testing.assert_array_almost_equal(_MM[i][1][1], MM_i[1][1])
        np.testing.assert_array_almost_equal(_MM[i][1][2], MM_i[1][2])
        np.testing.assert_array_almost_equal(_MM[i][2][0], MM_i[2][0])
        np.testing.assert_array_almost_equal(_MM[i][2][1], MM_i[2][1])
        np.testing.assert_array_almost_equal(_MM[i][2][2], MM_i[2][2])

        np.testing.assert_array_almost_equal(_iMM[i][0][0], iMM_i[0][0])
        np.testing.assert_array_almost_equal(_iMM[i][0][1], iMM_i[0][1])
        np.testing.assert_array_almost_equal(_iMM[i][0][2], iMM_i[0][2])
        np.testing.assert_array_almost_equal(_iMM[i][1][0], iMM_i[1][0])
        np.testing.assert_array_almost_equal(_iMM[i][1][1], iMM_i[1][1])
        np.testing.assert_array_almost_equal(_iMM[i][1][2], iMM_i[1][2])
        np.testing.assert_array_almost_equal(_iMM[i][2][0], iMM_i[2][0])
        np.testing.assert_array_almost_equal(_iMM[i][2][1], iMM_i[2][1])
        np.testing.assert_array_almost_equal(_iMM[i][2][2], iMM_i[2][2])

    return 1


def test_Mesh_NO4_elements_CT_QUAD():
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO4_elements_CT_QUAD} ...... ", flush=True)

    mesh_1 = MeshGenerator('crazy_periodic', c=0.25)([3, 2, 4], EDM='debug')
    mesh_2 = MeshGenerator('crazy_periodic')([2, 3, 4], EDM='debug')

    if rAnk == mAster_rank:
        ii, jj, kk = random.randint(1,5), random.randint(2,4), random.randint(2,3)
        quad_type = ['Gauss', 'Lobatto'][random.randint(0,1)]
    else:
        ii, jj, kk = None, None, None
        quad_type = None
    ii, jj, kk = cOmm.bcast([ii, jj, kk], root=mAster_rank)
    quad_type = cOmm.bcast(quad_type, root=mAster_rank)
    quad_degree = [ii, jj, kk]
    quad_nodes, quad_weights = Quadrature(quad_degree, category=quad_type).quad

    r, s, t = np.meshgrid(*quad_nodes, indexing='ij')

    for m in (mesh_1, mesh_2):
        _mapping = m.elements.coordinate_transformation.mapping(r, s, t)
        _X = m.elements.coordinate_transformation.X(r, s, t)
        _Y = m.elements.coordinate_transformation.Y(r, s, t)
        _Z = m.elements.coordinate_transformation.Z(r, s, t)
        _JM = m.elements.coordinate_transformation.Jacobian_matrix(r, s, t)
        _J00 = m.elements.coordinate_transformation.J00(r, s, t)
        _J01 = m.elements.coordinate_transformation.J01(r, s, t)
        _J02 = m.elements.coordinate_transformation.J02(r, s, t)
        _J10 = m.elements.coordinate_transformation.J10(r, s, t)
        _J11 = m.elements.coordinate_transformation.J11(r, s, t)
        _J12 = m.elements.coordinate_transformation.J12(r, s, t)
        _J20 = m.elements.coordinate_transformation.J20(r, s, t)
        _J21 = m.elements.coordinate_transformation.J21(r, s, t)
        _J22 = m.elements.coordinate_transformation.J22(r, s, t)
        _J = m.elements.coordinate_transformation.Jacobian(r, s, t, J=_JM)
        _J_ = m.elements.coordinate_transformation.Jacobian(r, s, t)
        _M = m.elements.coordinate_transformation.metric(r, s, t, detJ=_J)
        _M_ = m.elements.coordinate_transformation.metric(r, s, t)
        _MM = m.elements.coordinate_transformation.metric_matrix(r, s, t, J=_JM)
        _MM_ = m.elements.coordinate_transformation.metric_matrix(r, s, t)
        _iJM = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t, J=_JM)
        _iJM_ = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t)
        _iJ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t, iJ=_iJM)
        _iJ_ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t)
        _iMM = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t, iJ=_iJM)
        _iMM_ = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t)

        for i in m.elements:
            np.testing.assert_array_equal(_J[i], _J_[i])
            np.testing.assert_array_equal(_M[i], _M_[i])
            np.testing.assert_array_equal(_MM[i], _MM_[i])
            np.testing.assert_array_equal(_iJM[i], _iJM_[i])
            np.testing.assert_array_equal(_iJ[i], _iJ_[i])
            np.testing.assert_array_equal(_iMM[i], _iMM_[i])

        Q3_mapping = m.elements.coordinate_transformation.QUAD_3d.mapping(quad_degree, quad_type)
        Q3_X = m.elements.coordinate_transformation.QUAD_3d.X(quad_degree, quad_type)
        Q3_Y = m.elements.coordinate_transformation.QUAD_3d.Y(quad_degree, quad_type)
        Q3_Z = m.elements.coordinate_transformation.QUAD_3d.Z(quad_degree, quad_type)
        Q3_JM = m.elements.coordinate_transformation.QUAD_3d.Jacobian_matrix(quad_degree, quad_type)
        Q3_J00 = m.elements.coordinate_transformation.QUAD_3d.J00(quad_degree, quad_type)
        Q3_J01 = m.elements.coordinate_transformation.QUAD_3d.J01(quad_degree, quad_type)
        Q3_J02 = m.elements.coordinate_transformation.QUAD_3d.J02(quad_degree, quad_type)
        Q3_J10 = m.elements.coordinate_transformation.QUAD_3d.J10(quad_degree, quad_type)
        Q3_J11 = m.elements.coordinate_transformation.QUAD_3d.J11(quad_degree, quad_type)
        Q3_J12 = m.elements.coordinate_transformation.QUAD_3d.J12(quad_degree, quad_type)
        Q3_J20 = m.elements.coordinate_transformation.QUAD_3d.J20(quad_degree, quad_type)
        Q3_J21 = m.elements.coordinate_transformation.QUAD_3d.J21(quad_degree, quad_type)
        Q3_J22 = m.elements.coordinate_transformation.QUAD_3d.J22(quad_degree, quad_type)
        Q3_J = m.elements.coordinate_transformation.QUAD_3d.Jacobian(quad_degree, quad_type)
        Q3_M = m.elements.coordinate_transformation.QUAD_3d.metric(quad_degree, quad_type)
        Q3_MM = m.elements.coordinate_transformation.QUAD_3d.metric_matrix(quad_degree, quad_type)
        Q3_iJM = m.elements.coordinate_transformation.QUAD_3d.inverse_Jacobian_matrix(quad_degree, quad_type)
        Q3_iJ = m.elements.coordinate_transformation.QUAD_3d.inverse_Jacobian(quad_degree, quad_type)
        Q3_iMM = m.elements.coordinate_transformation.QUAD_3d.inverse_metric_matrix(quad_degree, quad_type)

        for i in m.elements:
            np.testing.assert_array_almost_equal(_mapping[i], Q3_mapping[i])
            np.testing.assert_array_almost_equal(_X[i], Q3_X[i])
            np.testing.assert_array_almost_equal(_Y[i], Q3_Y[i])
            np.testing.assert_array_almost_equal(_Z[i], Q3_Z[i])
            for j in range(3):
                for k in range(3):
                    np.testing.assert_array_almost_equal(_JM[i][j][k], Q3_JM[i][j][k])
            np.testing.assert_array_almost_equal(_J00[i], Q3_J00[i])
            np.testing.assert_array_almost_equal(_J01[i], Q3_J01[i])
            np.testing.assert_array_almost_equal(_J02[i], Q3_J02[i])
            np.testing.assert_array_almost_equal(_J10[i], Q3_J10[i])
            np.testing.assert_array_almost_equal(_J11[i], Q3_J11[i])
            np.testing.assert_array_almost_equal(_J12[i], Q3_J12[i])
            np.testing.assert_array_almost_equal(_J20[i], Q3_J20[i])
            np.testing.assert_array_almost_equal(_J21[i], Q3_J21[i])
            np.testing.assert_array_almost_equal(_J22[i], Q3_J22[i])
            np.testing.assert_array_almost_equal(_J[i], Q3_J[i])
            np.testing.assert_array_almost_equal(_M[i], Q3_M[i])
            np.testing.assert_array_almost_equal(_MM[i], Q3_MM[i])
            np.testing.assert_array_almost_equal(_iJM[i], Q3_iJM[i])
            np.testing.assert_array_almost_equal(_iJ[i], Q3_iJ[i])
            np.testing.assert_array_almost_equal(_iMM[i], Q3_iMM[i])

    r = r.ravel('F')
    s = s.ravel('F')
    t = t.ravel('F')
    for m in (mesh_1, mesh_2):
        _mapping = m.elements.coordinate_transformation.mapping(r, s, t)
        _X = m.elements.coordinate_transformation.X(r, s, t)
        _Y = m.elements.coordinate_transformation.Y(r, s, t)
        _Z = m.elements.coordinate_transformation.Z(r, s, t)
        _JM = m.elements.coordinate_transformation.Jacobian_matrix(r, s, t)
        _J00 = m.elements.coordinate_transformation.J00(r, s, t)
        _J01 = m.elements.coordinate_transformation.J01(r, s, t)
        _J02 = m.elements.coordinate_transformation.J02(r, s, t)
        _J10 = m.elements.coordinate_transformation.J10(r, s, t)
        _J11 = m.elements.coordinate_transformation.J11(r, s, t)
        _J12 = m.elements.coordinate_transformation.J12(r, s, t)
        _J20 = m.elements.coordinate_transformation.J20(r, s, t)
        _J21 = m.elements.coordinate_transformation.J21(r, s, t)
        _J22 = m.elements.coordinate_transformation.J22(r, s, t)
        _J = m.elements.coordinate_transformation.Jacobian(r, s, t, J=_JM)
        _J_ = m.elements.coordinate_transformation.Jacobian(r, s, t)
        _M = m.elements.coordinate_transformation.metric(r, s, t, detJ=_J)
        _M_ = m.elements.coordinate_transformation.metric(r, s, t)
        _MM = m.elements.coordinate_transformation.metric_matrix(r, s, t, J=_JM)
        _MM_ = m.elements.coordinate_transformation.metric_matrix(r, s, t)
        _iJM = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t, J=_JM)
        _iJM_ = m.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, t)
        _iJ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t, iJ=_iJM)
        _iJ_ = m.elements.coordinate_transformation.inverse_Jacobian(r, s, t)
        _iMM = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t, iJ=_iJM)
        _iMM_ = m.elements.coordinate_transformation.inverse_metric_matrix(r, s, t)

        for i in m.elements:
            np.testing.assert_array_equal(_J[i], _J_[i])
            np.testing.assert_array_equal(_M[i], _M_[i])
            np.testing.assert_array_equal(_MM[i], _MM_[i])
            np.testing.assert_array_equal(_iJM[i], _iJM_[i])
            np.testing.assert_array_equal(_iJ[i], _iJ_[i])
            np.testing.assert_array_equal(_iMM[i], _iMM_[i])

        Q3_mapping = m.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, quad_type)
        Q3_X = m.elements.coordinate_transformation.QUAD_1d.X(quad_degree, quad_type)
        Q3_Y = m.elements.coordinate_transformation.QUAD_1d.Y(quad_degree, quad_type)
        Q3_Z = m.elements.coordinate_transformation.QUAD_1d.Z(quad_degree, quad_type)
        Q3_JM = m.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, quad_type)
        Q3_J00 = m.elements.coordinate_transformation.QUAD_1d.J00(quad_degree, quad_type)
        Q3_J01 = m.elements.coordinate_transformation.QUAD_1d.J01(quad_degree, quad_type)
        Q3_J02 = m.elements.coordinate_transformation.QUAD_1d.J02(quad_degree, quad_type)
        Q3_J10 = m.elements.coordinate_transformation.QUAD_1d.J10(quad_degree, quad_type)
        Q3_J11 = m.elements.coordinate_transformation.QUAD_1d.J11(quad_degree, quad_type)
        Q3_J12 = m.elements.coordinate_transformation.QUAD_1d.J12(quad_degree, quad_type)
        Q3_J20 = m.elements.coordinate_transformation.QUAD_1d.J20(quad_degree, quad_type)
        Q3_J21 = m.elements.coordinate_transformation.QUAD_1d.J21(quad_degree, quad_type)
        Q3_J22 = m.elements.coordinate_transformation.QUAD_1d.J22(quad_degree, quad_type)
        Q3_J = m.elements.coordinate_transformation.QUAD_1d.Jacobian(quad_degree, quad_type)
        Q3_M = m.elements.coordinate_transformation.QUAD_1d.metric(quad_degree, quad_type)
        Q3_MM = m.elements.coordinate_transformation.QUAD_1d.metric_matrix(quad_degree, quad_type)
        Q3_iJM = m.elements.coordinate_transformation.QUAD_1d.inverse_Jacobian_matrix(quad_degree, quad_type)
        Q3_iJ = m.elements.coordinate_transformation.QUAD_1d.inverse_Jacobian(quad_degree, quad_type)
        Q3_iMM = m.elements.coordinate_transformation.QUAD_1d.inverse_metric_matrix(quad_degree, quad_type)

        for i in m.elements:
            np.testing.assert_array_almost_equal(_mapping[i], Q3_mapping[i])
            np.testing.assert_array_almost_equal(_X[i], Q3_X[i])
            np.testing.assert_array_almost_equal(_Y[i], Q3_Y[i])
            np.testing.assert_array_almost_equal(_Z[i], Q3_Z[i])
            for j in range(3):
                for k in range(3):
                    np.testing.assert_array_almost_equal(_JM[i][j][k], Q3_JM[i][j][k])
            np.testing.assert_array_almost_equal(_J00[i], Q3_J00[i])
            np.testing.assert_array_almost_equal(_J01[i], Q3_J01[i])
            np.testing.assert_array_almost_equal(_J02[i], Q3_J02[i])
            np.testing.assert_array_almost_equal(_J10[i], Q3_J10[i])
            np.testing.assert_array_almost_equal(_J11[i], Q3_J11[i])
            np.testing.assert_array_almost_equal(_J12[i], Q3_J12[i])
            np.testing.assert_array_almost_equal(_J20[i], Q3_J20[i])
            np.testing.assert_array_almost_equal(_J21[i], Q3_J21[i])
            np.testing.assert_array_almost_equal(_J22[i], Q3_J22[i])
            np.testing.assert_array_almost_equal(_J[i], Q3_J[i])
            np.testing.assert_array_almost_equal(_M[i], Q3_M[i])
            np.testing.assert_array_almost_equal(_MM[i], Q3_MM[i])
            np.testing.assert_array_almost_equal(_iJM[i], Q3_iJM[i])
            np.testing.assert_array_almost_equal(_iJ[i], Q3_iJ[i])
            np.testing.assert_array_almost_equal(_iMM[i], Q3_iMM[i])

    return 1


def test_Mesh_NO5_mesh_trace_topology():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO5_mesh_trace_topology} ...... ", flush=True)

    MID = list(DomainInputAllocator.___defined_DI___().keys())
    if rAnk == mAster_rank:
        __ = random.sample(range(0,len(MID)), 2)
        meshes = [MID[i] for i in __]
        II = random.randint(3,4) # [II, JJ, KK] element layout
        JJ = random.randint(2,5) # [II, JJ, KK] element layout
        KK = random.randint(1,4) # [II, JJ, KK] element layout
    else:
        meshes = None
        II, JJ, KK = None, None, None
    II, JJ, KK = cOmm.bcast([II, JJ, KK], root=mAster_rank)
    meshes = cOmm.bcast(meshes, root=mAster_rank)

    for mid in meshes:
        # ... generate meshes ...
        if mid in ('crazy', 'crazy_periodic'):
            if rAnk == mAster_rank:
                c = random.uniform(0, 0.3)
            else:
                c = None
            c = cOmm.bcast(c, root=mAster_rank)
            mesh = MeshGenerator(mid, c=c)([II, JJ, KK], EDM='debug')
        else:
            try:
                mesh = MeshGenerator(mid)([II, JJ, KK], EDM='debug')
            except ThreeDimensionalTransfiniteInterpolationError:
                mesh = MeshGenerator('crazy')([II, JJ, KK], EDM='debug')

        elements = mesh.elements
        SD = list()

        MAP = mesh.trace.elements.map
        for ele_i in MAP:
            for i in MAP[ele_i]:
                assert i in mesh.trace.elements

        for i in mesh.trace.elements:
            e = mesh.trace.elements[i]
            assert e.i == i
            shared_with_core = e.shared_with_core
            assert e.CHARACTERISTIC_element in elements
            if shared_with_core is None:
                pass
            else:
                SD.extend([rAnk, shared_with_core])

            if e.IS.on_mesh_boundary:
                assert e.positions[1] in mesh.domain.boundaries.names
            if e.IS.on_periodic_boundary:
                assert not e.IS.on_mesh_boundary
                assert e.positions[1][0] in '0123456789'

        SD = cOmm.gather(SD, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            sd = list()
            for SDi in SD:
                sd.extend(SDi)
            sd_SET =set(sd)
            for i in sd_SET:
                assert sd.count(i) % 2 == 0

    return 1


def test_Mesh_NO5a_mesh_trace_CT():
    """
    Unittests for the mesh - trace elements - CT.
    """
    if rAnk == mAster_rank:
        print("ttt {test_Mesh_NO5a_mesh_trace_CT} ...... ", flush=True)


    if rAnk == mAster_rank:
        el1 = random.randint(1,4)
        el2 = random.randint(1,3)
        el3 = random.randint(2,3)
        c = random.uniform(0., 0.25)
        if c < 0.1:
            c = 0
    else:
        el1, el2, el3, c = None, None, None, None
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)
    M1 = MeshGenerator('crazy_periodic', c=c)([el1, el2, el3])


    if rAnk == mAster_rank:
        el1 = random.randint(1,4)
        el2 = random.randint(1,3)
        el3 = random.randint(2,3)
        c = random.uniform(0., 0.25)
        if c < 0.1:
            c = 0
    else:
        el1, el2, el3, c = None, None, None, None
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)
    M2 = MeshGenerator('crazy', c=c)([el1, el2, el3])

    for M in (M1, M2):

        tes = M.trace.elements

        xi = np.random.rand(3,3)
        et = np.random.rand(3,3)
        sg = np.random.rand(3,3)

        JM = tes.coordinate_transformation.Jacobian_matrix(xi, et, sg)
        iJM = tes.coordinate_transformation.inverse_Jacobian_matrix(xi, et, sg)
        MM = tes.coordinate_transformation.metric_matrix(xi, et, sg)
        MT = tes.coordinate_transformation.metric(xi, et, sg)
        UNV = tes.coordinate_transformation.unit_normal_vector(xi, et, sg)

        for i in tes:
            te = tes[i]
            side = te.CHARACTERISTIC_side

            if side in 'NS':
                _xi_eta_sigma_ = [et, sg]
            elif side in 'WE':
                _xi_eta_sigma_ = [xi, sg]
            elif side in 'BF':
                _xi_eta_sigma_ = [xi, et]
            else:
                raise Exception()

            jm = te.coordinate_transformation.Jacobian_matrix(*_xi_eta_sigma_)
            ijm = te.coordinate_transformation.inverse_Jacobian_matrix(*_xi_eta_sigma_)
            mm = te.coordinate_transformation.metric_matrix(*_xi_eta_sigma_)
            mt = te.coordinate_transformation.metric(*_xi_eta_sigma_)
            unv = te.coordinate_transformation.unit_normal_vector(*_xi_eta_sigma_)

            np.testing.assert_almost_equal(MT[i], mt)

            MMi = MM[i]
            np.testing.assert_almost_equal(MMi[0][0], mm[0][0])
            np.testing.assert_almost_equal(MMi[0][1], mm[0][1])
            np.testing.assert_almost_equal(MMi[1][0], mm[1][0])
            np.testing.assert_almost_equal(MMi[1][1], mm[1][1])

            JMi = JM[i]
            np.testing.assert_almost_equal(JMi[0][0], jm[0][0])
            np.testing.assert_almost_equal(JMi[0][1], jm[0][1])
            np.testing.assert_almost_equal(JMi[1][0], jm[1][0])
            np.testing.assert_almost_equal(JMi[1][1], jm[1][1])
            np.testing.assert_almost_equal(JMi[2][0], jm[2][0])
            np.testing.assert_almost_equal(JMi[2][1], jm[2][1])

            iJMi = iJM[i]
            np.testing.assert_almost_equal(iJMi[0][0], ijm[0][0])
            np.testing.assert_almost_equal(iJMi[0][1], ijm[0][1])
            np.testing.assert_almost_equal(iJMi[0][2], ijm[0][2])
            np.testing.assert_almost_equal(iJMi[1][0], ijm[1][0])
            np.testing.assert_almost_equal(iJMi[1][1], ijm[1][1])
            np.testing.assert_almost_equal(iJMi[1][2], ijm[1][2])

            UNVi = UNV[i]
            np.testing.assert_almost_equal(UNVi[0], unv[0])
            np.testing.assert_almost_equal(UNVi[1], unv[1])
            np.testing.assert_almost_equal(UNVi[2], unv[2])



    return 1

def test_Mesh_NO6_transfinite():
    """Unittests for the mesh."""
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO6_transfinite} ...... ", flush=True)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    try:

        mesh = MeshGenerator('psc')([4,2,2])

        space = SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',4), ('Lobatto',4)])
        FC = FormCaller(mesh, space)

        scalar = FC('scalar', p)
        vector = FC('vector', (u,v,w))

        f0 = FC('0-f', is_hybrid=False)
        f1 = FC('1-f', is_hybrid=False)
        f2 = FC('2-f', is_hybrid=False)
        f3 = FC('3-f', is_hybrid=False)

        f0.TW.func.body = scalar
        f0.TW.___DO_push_all_to_instant___(0)
        f0.discretize()
        assert f0.error.L() < 0.0022
        f1.TW.func.body = vector
        f1.TW.___DO_push_all_to_instant___(0)
        f1.discretize()
        assert f1.error.L() < 0.0043
        f2.TW.func.body = vector
        f2.TW.___DO_push_all_to_instant___(0)
        f2.discretize()
        assert f2.error.L() < 0.0048
        f3.TW.func.body = scalar
        f3.TW.___DO_push_all_to_instant___(0)
        f3.discretize()
        assert f3.error.L() < 0.003

    except ThreeDimensionalTransfiniteInterpolationError:

        if rAnk == mAster_rank:
            print("   ~ Transfinite test SKIPPED.", flush=True)

    return 1

def test_Mesh_NO7_boundaries():
    """Unittests for the mesh."""
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO7_boundaries} ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic')([3, 3, 3], EDM=None, show_info=False)

    DB = mesh.domain.boundaries
    MB = mesh.boundaries

    DBN = DB.names
    MBN = MB.names

    # below, we test that at domain.boundaries, the periodic boundaries are included while in mesh.boundaries they are not.
    assert DBN == ('North', 'South', 'West', 'East', 'Back', 'Front')
    assert MBN == tuple()

    return 1

def test_Mesh_NO8_Mesh_SubGeometry_perpendicular_slice_object():
    """Unittests for the mesh.

    Also used to show how to generate Mesh_SubGeometry.
    """
    if rAnk == mAster_rank:
        print(">>> {test_Mesh_NO8_Mesh_SubGeometry_perpendicular_slice_object} ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic')([3, 3, 3], EDM=None, show_info=False)
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)


    R = mesh.domain.regions['R:R']
    RSG = R.sub_geometry
    RSG_PSO = RSG.make_a_perpendicular_slice_object_on(r=0.5)
    MSG_PSO = mesh.sub_geometry.make_a_perpendicular_slice_object_on(RSG_PSO)


    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return x + np.sin(2*np.pi*y)*np.sin(2*np.pi*z) + t/2
    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f2 = FC('2-f', is_hybrid=False)
    f3 = FC('3-f', is_hybrid=False)
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()

    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()

    f2.TW.func.body = vector
    f2.TW.do.push_all_to_instant(0)
    f2.discretize()

    f3.TW.func.body = scalar
    f3.TW.do.push_all_to_instant(0)
    f3.discretize()

    f0.visualize.matplot.perpendicular_slice(MSG_PSO, usetex=False, saveto='No8_perpendicular_slice_object_f0.pdf')
    f1.visualize.matplot.perpendicular_slice(MSG_PSO, usetex=False, saveto='No8_perpendicular_slice_object_f1.pdf')
    f2.visualize.matplot.perpendicular_slice(MSG_PSO, usetex=False, saveto='No8_perpendicular_slice_object_f2.pdf')
    f3.visualize.matplot.perpendicular_slice(MSG_PSO, usetex=False, saveto='No8_perpendicular_slice_object_f3.pdf')

    if rAnk == mAster_rank:
        os.remove("No8_perpendicular_slice_object_f0.pdf")
        os.remove("No8_perpendicular_slice_object_f1_0th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f1_1th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f1_2th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f2_0th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f2_1th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f2_2th_component.pdf")
        os.remove("No8_perpendicular_slice_object_f3.pdf")

    return 1





def test_Mesh_NO9_edge_node_mesh():
    if rAnk == mAster_rank:
        print("ENM {test_Mesh_NO9_edge_node_mesh} ...... ", flush=True)

    if rAnk == mAster_rank:
        LOAD = random.randint(50, 1000)
    else:
        LOAD = None
    LOAD = cOmm.bcast(LOAD, root=mAster_rank)
    mesh = random_mesh_of_elements_around(LOAD, mesh_pool=['bridge_arch_cracked', ], EDM_pool=['chaotic', ])
    # add to mesh_pool to test it with more meshes.
    MN = mesh.node
    MNE = MN.elements
    ME = mesh.edge
    MEE = ME.elements

    # ---- topology test: the locations of node elements -------------------------------------------
    locations = MNE._locations_
    for node in locations:
        assert len(locations[node]) == len(set(locations[node])), f"a trivial check!"

        elements = list()
        for loc in locations[node]:
            if loc[0] in '0123456789':
                elements.append(loc[:-3])

        assert len(elements) == len(set(elements)), \
            f"a node element cannot be two corners of one mesh element unless it is a fully periodic domain" \
            f"of one 1 mesh element along an axis, which is not allowed!"

    if sIze > 1: # only need to do this check when use >1 cores.

        for i in range(sIze):
            LOCATIONS = cOmm.bcast(locations, root=i)

            if rAnk != i: # do the check
                for node in LOCATIONS:
                    if node in locations:
                        assert set(locations[node]) == set(LOCATIONS[node]), \
                            f"location[{node}] = {locations[node]} in core #{rAnk} is not equal to" \
                            f"location[{node}] = {LOCATIONS[node]} in core #{i}."

    # ---- topology test: the locations of edge elements -------------------------------------------
    locations = MEE._locations_
    for edge in locations:
        assert len(locations[edge]) == len(set(locations[edge])), f"a trivial check!"

        elements = list()
        for loc in locations[edge]:
            if loc[0] in '0123456789':
                elements.append(loc[:-2])

        assert len(elements) == len(set(elements)), \
            f"an edge element cannot be two corner edges of one mesh element unless it is a fully periodic domain" \
            f"of one 1 mesh element along an axis, which is not allowed!"

    if sIze > 1: # only need to do this check when use >1 cores.

        for i in range(sIze):
            LOCATIONS = cOmm.bcast(locations, root=i)

            if rAnk != i: # do the check
                for edge in LOCATIONS:
                    if edge in locations:
                        assert set(locations[edge]) == set(LOCATIONS[edge]), \
                            f"location[{edge}] = {locations[edge]} in core #{rAnk} is not equal to" \
                            f"location[{edge}] = {LOCATIONS[edge]} in core #{i}."

    return  1



if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\TESTS\unittest_mesh.py
    test_Mesh_NO9_edge_node_mesh()