# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from objects.CSCG._2d.master import MeshGenerator
from objects.CSCG._2d.mesh.domain.inputs.allocator import DomainInputFinder
import random
from screws.quadrature import Quadrature



def test_Mesh_NO1_mesh_topology():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("+++ [test_Mesh_NO1_mesh_topology] ...... ", flush=True)

    mesh = MeshGenerator('chp2')([3, 4], EDM='debug')
    MAP = mesh.elements.map
    if 0 in MAP: assert MAP[0] == ('Upper', 1, 'Left', 3)
    if 1 in MAP: assert MAP[1] == (0, 2, 'Left', 4)
    if 73 in MAP: assert MAP[73] == (72, 74, 'Internal', 76)
    if 33 in MAP: assert MAP[33] == (23, 34, 30, 48)
    np.testing.assert_array_equal(mesh.elements.spacing['R:R_UL'][1],
                                  np.array([0.  ,0.25, 0.5 , 0.75, 1.  ]))

    mesh = MeshGenerator('chp1')({'R:Ru':[4, 3],
                                  'R:Rl':[5, 3],
                                  'R:Rd':[[1,2], 3],
                                  'R:Rr':[[1,2,1], 3],}, EDM='debug')
    MAP = mesh.elements.map
    if 0 in MAP: assert MAP[0] == (31, 1, 'Internal', 4)
    if 12 in MAP: assert MAP[12] == (3, 13, 'Internal', 15)
    if 26 in MAP: assert MAP[26] == (25, 37, 24, 'Down')
    if 32 in MAP: assert MAP[32] == (24, 33, 27, 37)
    if 39 in MAP: assert MAP[39] == (38, 40, 34, 'Left')
    if 19 in MAP: assert MAP[19] == (18, 20, 16, 'Right')

    mesh = MeshGenerator('chp2')([3, 4], EDM='debug')
    for rn in mesh.domain.regions.names:
        R = mesh.domain.regions[rn]
        if rn in ['R:R_DR', 'R:R_UL', 'R:R_DL', 'R:R_UR']:
            assert R.type_wrt_metric.mark == 'orthogonal:UD0.64644661_LR0.64644661'
        else:
            assert R.type_wrt_metric.mark[:8] == 'chaotic:'


    mesh = MeshGenerator('cic')([3, 4], EDM='debug')
    for rn in mesh.domain.regions.names:
        R = mesh.domain.regions[rn]
        if rn == 'R:Ri':
            assert R.type_wrt_metric.mark == \
                   'parallelogram:angleL1.57079633_lenL1.50000000_A1.57079633A_lenU0.75000000'
        elif rn == 'R:Ro':
            assert R.type_wrt_metric.mark == \
                   'parallelogram:angleL4.71238898_lenL1.50000000_A1.57079633A_lenU2.25000000'
        else:
            assert R.type_wrt_metric.mark[:8] == 'chaotic:'

    mesh = MeshGenerator('crazy_periodic', c=0.3)([3, 4], EDM='debug')
    MAP = mesh.elements.map
    if 0 in MAP: assert MAP[0] == (2, 1, 9, 3)
    if 1 in MAP: assert MAP[1] == (0, 2, 10, 4)
    if 2 in MAP: assert MAP[2] == (1, 0, 11, 5)
    if 3 in MAP: assert MAP[3] == (5, 4, 0, 6)
    if 4 in MAP: assert MAP[4] == (3, 5, 1, 7)
    if 5 in MAP: assert MAP[5] == (4, 3, 2, 8)
    if 7 in MAP: assert MAP[7] == (6, 8, 4, 10)
    if 9 in MAP: assert MAP[9] == (11, 10, 6, 0)
    if 10 in MAP: assert MAP[10] == (9, 11, 7, 1)
    if 11 in MAP: assert MAP[11] == (10, 9, 8, 2)

    mesh = MeshGenerator('quadrangle')([3, 4], EDM=None)
    for i in mesh.elements:
        element = mesh.elements[i]
        mark = element.type_wrt_metric.mark
        assert mark[:13] == 'Parallelogram', "error!"
    mesh = MeshGenerator('quadrangle', p_UL=(0,0), p_DL=(1,0), p_UR=(0,1), p_DR=(1,1))(
        [3, 4], EDM=None)
    for i in mesh.elements:
        element = mesh.elements[i]
        mark = element.type_wrt_metric.mark
        assert mark[:4] == 'Orth', "error!"
    mesh = MeshGenerator('quadrangle', p_UL=(1,0), p_DL=(2,1), p_UR=(0,1), p_DR=(1,2))(
        [3, 4], EDM=None)
    for i in mesh.elements:
        element = mesh.elements[i]
        mark = element.type_wrt_metric.mark
        assert mark[:13] == 'Parallelogram', "error!"
    
    return 1


def test_Mesh_NO2_mesh_coordinate_transformation():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("+++ [test_Mesh_NO2_mesh_coordinate_transformation] ...... ", flush=True)


    MID = list(DomainInputFinder.___defined_DI___().keys())

    if rAnk == mAster_rank:
        __ = random.sample(range(0,len(MID)), 4)
        meshes = [MID[i] for i in __]
        II = random.randint(3,4) # [II, JJ] element layout
        JJ = random.randint(2,5) # [II, JJ] element layout
    else:
        meshes = None
        II, JJ = None, None
    II, JJ = cOmm.bcast([II, JJ], root=mAster_rank)
    meshes = cOmm.bcast(meshes, root=mAster_rank)
    for mid in meshes:
        # ... generate meshes ...
        if mid in ('crazy', 'crazy_periodic'):
            if rAnk == mAster_rank:
                c = random.uniform(0, 0.3)
            else:
                c = None
            c = cOmm.bcast(c, root=mAster_rank)
            mesh = MeshGenerator(mid, c=c)([II, JJ], EDM='debug')
        else:
            mesh = MeshGenerator(mid)([II, JJ], EDM='debug')
        # ... generate r, s...
        if rAnk == mAster_rank:
            r = np.linspace(-0.95, 0.975, random.randint(2,8))
            s = np.linspace(random.uniform(-0.99, -0.9), random.uniform(0.85, 0.99), random.randint(1,7))
        else:
            r, s = None, None
        r, s = cOmm.bcast([r, s], root=mAster_rank)
        #... now lets check the coordinate transformation ...

        r, s = np.meshgrid(r, s, indexing='ij')

        mesh.___TEST_MODE___ = True
        mesh.___DEPRECATED_ct___.evaluated_at(r, s)

        mapping = mesh.___DEPRECATED_ct___.mapping
        JM = mesh.___DEPRECATED_ct___.Jacobian_matrix
        J = mesh.___DEPRECATED_ct___.Jacobian
        iJM = mesh.___DEPRECATED_ct___.inverse_Jacobian_matrix
        iJ = mesh.___DEPRECATED_ct___.inverse_Jacobian

        M = mesh.___DEPRECATED_ct___.metric
        MM = mesh.___DEPRECATED_ct___.metric_matrix
        iMM = mesh.___DEPRECATED_ct___.inverse_metric_matrix

        _mapping = mesh.elements.coordinate_transformation.mapping(r, s)
        _X = mesh.elements.coordinate_transformation.X(r, s)
        _Y = mesh.elements.coordinate_transformation.Y(r, s)
        _JM = mesh.elements.coordinate_transformation.Jacobian_matrix(r, s)
        _J00 = mesh.elements.coordinate_transformation.J00(r, s)
        _J01 = mesh.elements.coordinate_transformation.J01(r, s)
        _J10 = mesh.elements.coordinate_transformation.J10(r, s)
        _J11 = mesh.elements.coordinate_transformation.J11(r, s)
        _J = mesh.elements.coordinate_transformation.Jacobian(r, s, J=_JM)
        _M = mesh.elements.coordinate_transformation.metric(r, s, detJ=_J)
        _MM = mesh.elements.coordinate_transformation.metric_matrix(r, s, J=_JM)
        _iJM = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, J=_JM)
        _iJ = mesh.elements.coordinate_transformation.inverse_Jacobian(r, s, iJ=_iJM)
        _iMM = mesh.elements.coordinate_transformation.inverse_metric_matrix(r, s, iJ=_iJM)

        for i in mesh.elements.indices:
            ei = mesh.elements[i]

            mapping_i = ei.coordinate_transformation.mapping(r, s)
            X = ei.coordinate_transformation.X(r, s)
            Y = ei.coordinate_transformation.Y(r, s)

            np.testing.assert_array_almost_equal(mapping[0][i], X)
            np.testing.assert_array_almost_equal(mapping[1][i], Y)
            np.testing.assert_array_almost_equal(_mapping[i][0], X)
            np.testing.assert_array_almost_equal(_mapping[i][1], Y)
            np.testing.assert_array_almost_equal(_X[i], X)
            np.testing.assert_array_almost_equal(_Y[i], Y)

            np.testing.assert_array_almost_equal(mapping[0][i], mapping_i[0])
            np.testing.assert_array_almost_equal(mapping[1][i], mapping_i[1])

            JM_i = ei.coordinate_transformation.Jacobian_matrix(r, s)

            np.testing.assert_array_almost_equal(JM[0][0][i], JM_i[0][0])
            np.testing.assert_array_almost_equal(JM[0][1][i], JM_i[0][1])
            np.testing.assert_array_almost_equal(JM[1][0][i], JM_i[1][0])
            np.testing.assert_array_almost_equal(JM[1][1][i], JM_i[1][1])
            np.testing.assert_array_almost_equal(_JM[i][0][0], JM_i[0][0])
            np.testing.assert_array_almost_equal(_JM[i][0][1], JM_i[0][1])
            np.testing.assert_array_almost_equal(_JM[i][1][0], JM_i[1][0])
            np.testing.assert_array_almost_equal(_JM[i][1][1], JM_i[1][1])

            J00 = ei.coordinate_transformation.J00(r, s)
            J01 = ei.coordinate_transformation.J01(r, s)
            J10 = ei.coordinate_transformation.J10(r, s)
            J11 = ei.coordinate_transformation.J11(r, s)

            np.testing.assert_array_almost_equal(JM[0][0][i], J00)
            np.testing.assert_array_almost_equal(JM[0][1][i], J01)
            np.testing.assert_array_almost_equal(JM[1][0][i], J10)
            np.testing.assert_array_almost_equal(JM[1][1][i], J11)
            np.testing.assert_array_almost_equal(_J00[i], J00)
            np.testing.assert_array_almost_equal(_J01[i], J01)
            np.testing.assert_array_almost_equal(_J10[i], J10)
            np.testing.assert_array_almost_equal(_J11[i], J11)

            J0 = ei.coordinate_transformation.J0_(r, s)
            J1 = ei.coordinate_transformation.J1_(r, s)
            np.testing.assert_array_almost_equal(J0[0], J00)
            np.testing.assert_array_almost_equal(J0[1], J01)
            np.testing.assert_array_almost_equal(J1[0], J10)
            np.testing.assert_array_almost_equal(J1[1], J11)

            J_i = ei.coordinate_transformation.Jacobian(r, s)
            iJ_i = ei.coordinate_transformation.inverse_Jacobian(r, s)
            M_i = ei.coordinate_transformation.metric(r, s)

            np.testing.assert_array_almost_equal(J[i], J_i)
            np.testing.assert_array_almost_equal(_J[i], J_i)
            np.testing.assert_array_almost_equal(iJ[i], iJ_i)
            np.testing.assert_array_almost_equal(_iJ[i], iJ_i)
            np.testing.assert_array_almost_equal(M[i], M_i)
            np.testing.assert_array_almost_equal(_M[i], M_i)

            iJM_i = ei.coordinate_transformation.inverse_Jacobian_matrix(r, s)
            np.testing.assert_array_almost_equal(iJM[0][0][i], iJM_i[0][0])
            np.testing.assert_array_almost_equal(iJM[0][1][i], iJM_i[0][1])
            np.testing.assert_array_almost_equal(iJM[1][0][i], iJM_i[1][0])
            np.testing.assert_array_almost_equal(iJM[1][1][i], iJM_i[1][1])
            np.testing.assert_array_almost_equal(_iJM[i][0][0], iJM_i[0][0])
            np.testing.assert_array_almost_equal(_iJM[i][0][1], iJM_i[0][1])
            np.testing.assert_array_almost_equal(_iJM[i][1][0], iJM_i[1][0])
            np.testing.assert_array_almost_equal(_iJM[i][1][1], iJM_i[1][1])

            MM_i = ei.coordinate_transformation.metric_matrix(r, s)
            iMM_i = ei.coordinate_transformation.inverse_metric_matrix(r, s)
            np.testing.assert_array_almost_equal(MM[0][0][i], MM_i[0][0])
            np.testing.assert_array_almost_equal(MM[0][1][i], MM_i[0][1])
            np.testing.assert_array_almost_equal(MM[1][0][i], MM_i[1][0])
            np.testing.assert_array_almost_equal(MM[1][1][i], MM_i[1][1])
            np.testing.assert_array_almost_equal(_MM[i][0][0], MM_i[0][0])
            np.testing.assert_array_almost_equal(_MM[i][0][1], MM_i[0][1])
            np.testing.assert_array_almost_equal(_MM[i][1][0], MM_i[1][0])
            np.testing.assert_array_almost_equal(_MM[i][1][1], MM_i[1][1])

            np.testing.assert_array_almost_equal(iMM[0][0][i], iMM_i[0][0])
            np.testing.assert_array_almost_equal(iMM[0][1][i], iMM_i[0][1])
            np.testing.assert_array_almost_equal(iMM[1][0][i], iMM_i[1][0])
            np.testing.assert_array_almost_equal(iMM[1][1][i], iMM_i[1][1])
            np.testing.assert_array_almost_equal(_iMM[i][0][0], iMM_i[0][0])
            np.testing.assert_array_almost_equal(_iMM[i][0][1], iMM_i[0][1])
            np.testing.assert_array_almost_equal(_iMM[i][1][0], iMM_i[1][0])
            np.testing.assert_array_almost_equal(_iMM[i][1][1], iMM_i[1][1])

    return 1


def test_Mesh_NO3_mesh_coordinate_transformation_QUAD():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("+++ [test_Mesh_NO3_mesh_coordinate_transformation_QUAD] ...... ", flush=True)

    MID = list(MeshGenerator.___domain_input_statistic___().keys())

    if rAnk == mAster_rank:
        __ = random.sample(range(0,len(MID)), 3)
        meshes = [MID[i] for i in __]
        II = random.randint(3,4) # [II, JJ] element layout
        JJ = random.randint(2,3) # [II, JJ] element layout
    else:
        meshes = None
        II, JJ = None, None
    II, JJ = cOmm.bcast([II, JJ], root=mAster_rank)
    meshes = cOmm.bcast(meshes, root=mAster_rank)
    for mid in meshes:
        # ... generate meshes ...
        if mid in ('crazy', 'crazy_periodic'):
            if rAnk == mAster_rank:
                c = random.uniform(0, 0.3)
            else:
                c = None
            c = cOmm.bcast(c, root=mAster_rank)
            mesh = MeshGenerator(mid, c=c)([II, JJ], EDM='debug')
        else:
            mesh = MeshGenerator(mid)([II, JJ], EDM='debug')

        if rAnk == mAster_rank:
            quad_degree = [random.randint(3,5),random.randint(2,3)]
            quad_type = ['Gauss', 'Lobatto'][random.randint(0,1)]
        else:
            quad_degree, quad_type = None, None
        quad_degree, quad_type = cOmm.bcast([quad_degree, quad_type], root=mAster_rank)
        quad_nodes, quad_weights = Quadrature(quad_degree, category=quad_type).quad
        r, s = np.meshgrid(*quad_nodes, indexing='ij')

        _mapping = mesh.elements.coordinate_transformation.mapping(r, s)
        _X = mesh.elements.coordinate_transformation.X(r, s)
        _Y = mesh.elements.coordinate_transformation.Y(r, s)
        _JM = mesh.elements.coordinate_transformation.Jacobian_matrix(r, s)
        _J00 = mesh.elements.coordinate_transformation.J00(r, s)
        _J01 = mesh.elements.coordinate_transformation.J01(r, s)
        _J10 = mesh.elements.coordinate_transformation.J10(r, s)
        _J11 = mesh.elements.coordinate_transformation.J11(r, s)
        _J = mesh.elements.coordinate_transformation.Jacobian(r, s, J=_JM)
        _J_ = mesh.elements.coordinate_transformation.Jacobian(r, s)
        _M = mesh.elements.coordinate_transformation.metric(r, s, detJ=_J)
        _M_ = mesh.elements.coordinate_transformation.metric(r, s)
        _MM = mesh.elements.coordinate_transformation.metric_matrix(r, s, J=_JM)
        _MM_ = mesh.elements.coordinate_transformation.metric_matrix(r, s)
        _iJM = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, J=_JM)
        _iJM_ = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s)
        _iJ = mesh.elements.coordinate_transformation.inverse_Jacobian(r, s, iJ=_iJM)
        _iJ_ = mesh.elements.coordinate_transformation.inverse_Jacobian(r, s)
        _iMM = mesh.elements.coordinate_transformation.inverse_metric_matrix(r, s, iJ=_iJM)
        _iMM_ = mesh.elements.coordinate_transformation.inverse_metric_matrix(r, s)

        for i in mesh.elements:
            np.testing.assert_array_equal(_J[i], _J_[i])
            np.testing.assert_array_equal(_M[i], _M_[i])
            np.testing.assert_array_equal(_MM[i], _MM_[i])
            np.testing.assert_array_equal(_iJM[i], _iJM_[i])
            np.testing.assert_array_equal(_iJ[i], _iJ_[i])
            np.testing.assert_array_equal(_iMM[i], _iMM_[i])

        Q3_mapping = mesh.elements.coordinate_transformation.QUAD_2d.mapping(quad_degree, quad_type)
        Q3_X = mesh.elements.coordinate_transformation.QUAD_2d.X(quad_degree, quad_type)
        Q3_Y = mesh.elements.coordinate_transformation.QUAD_2d.Y(quad_degree, quad_type)
        Q3_JM = mesh.elements.coordinate_transformation.QUAD_2d.Jacobian_matrix(quad_degree, quad_type)
        Q3_J00 = mesh.elements.coordinate_transformation.QUAD_2d.J00(quad_degree, quad_type)
        Q3_J01 = mesh.elements.coordinate_transformation.QUAD_2d.J01(quad_degree, quad_type)
        Q3_J10 = mesh.elements.coordinate_transformation.QUAD_2d.J10(quad_degree, quad_type)
        Q3_J11 = mesh.elements.coordinate_transformation.QUAD_2d.J11(quad_degree, quad_type)
        Q3_J = mesh.elements.coordinate_transformation.QUAD_2d.Jacobian(quad_degree, quad_type)
        Q3_M = mesh.elements.coordinate_transformation.QUAD_2d.metric(quad_degree, quad_type)
        Q3_MM = mesh.elements.coordinate_transformation.QUAD_2d.metric_matrix(quad_degree, quad_type)
        Q3_iJM = mesh.elements.coordinate_transformation.QUAD_2d.inverse_Jacobian_matrix(quad_degree, quad_type)
        Q3_iJ = mesh.elements.coordinate_transformation.QUAD_2d.inverse_Jacobian(quad_degree, quad_type)
        Q3_iMM = mesh.elements.coordinate_transformation.QUAD_2d.inverse_metric_matrix(quad_degree, quad_type)

        for i in mesh.elements:
            np.testing.assert_array_almost_equal(_mapping[i], Q3_mapping[i])
            np.testing.assert_array_almost_equal(_X[i], Q3_X[i])
            np.testing.assert_array_almost_equal(_Y[i], Q3_Y[i])
            for j in range(2):
                for k in range(2):
                    np.testing.assert_array_almost_equal(_JM[i][j][k], Q3_JM[i][j][k])
            np.testing.assert_array_almost_equal(_J00[i], Q3_J00[i])
            np.testing.assert_array_almost_equal(_J01[i], Q3_J01[i])
            np.testing.assert_array_almost_equal(_J10[i], Q3_J10[i])
            np.testing.assert_array_almost_equal(_J11[i], Q3_J11[i])
            np.testing.assert_array_almost_equal(_J[i], Q3_J[i])
            np.testing.assert_array_almost_equal(_M[i], Q3_M[i])
            np.testing.assert_array_almost_equal(_MM[i], Q3_MM[i])
            np.testing.assert_array_almost_equal(_iJM[i], Q3_iJM[i])
            np.testing.assert_array_almost_equal(_iJ[i], Q3_iJ[i])
            np.testing.assert_array_almost_equal(_iMM[i], Q3_iMM[i])

        r = r.ravel('F')
        s = s.ravel('F')

        _mapping = mesh.elements.coordinate_transformation.mapping(r, s)
        _X = mesh.elements.coordinate_transformation.X(r, s)
        _Y = mesh.elements.coordinate_transformation.Y(r, s)
        _JM = mesh.elements.coordinate_transformation.Jacobian_matrix(r, s)
        _J00 = mesh.elements.coordinate_transformation.J00(r, s)
        _J01 = mesh.elements.coordinate_transformation.J01(r, s)
        _J10 = mesh.elements.coordinate_transformation.J10(r, s)
        _J11 = mesh.elements.coordinate_transformation.J11(r, s)
        _J = mesh.elements.coordinate_transformation.Jacobian(r, s, J=_JM)
        _J_ = mesh.elements.coordinate_transformation.Jacobian(r, s)
        _M = mesh.elements.coordinate_transformation.metric(r, s, detJ=_J)
        _M_ = mesh.elements.coordinate_transformation.metric(r, s)
        _MM = mesh.elements.coordinate_transformation.metric_matrix(r, s, J=_JM)
        _MM_ = mesh.elements.coordinate_transformation.metric_matrix(r, s)
        _iJM = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s, J=_JM)
        _iJM_ = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(r, s)
        _iJ = mesh.elements.coordinate_transformation.inverse_Jacobian(r, s, iJ=_iJM)
        _iJ_ = mesh.elements.coordinate_transformation.inverse_Jacobian(r, s)
        _iMM = mesh.elements.coordinate_transformation.inverse_metric_matrix(r, s, iJ=_iJM)
        _iMM_ = mesh.elements.coordinate_transformation.inverse_metric_matrix(r, s)

        for i in mesh.elements:
            np.testing.assert_array_equal(_J[i], _J_[i])
            np.testing.assert_array_equal(_M[i], _M_[i])
            np.testing.assert_array_equal(_MM[i], _MM_[i])
            np.testing.assert_array_equal(_iJM[i], _iJM_[i])
            np.testing.assert_array_equal(_iJ[i], _iJ_[i])
            np.testing.assert_array_equal(_iMM[i], _iMM_[i])

        Q3_mapping = mesh.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, quad_type)
        Q3_X = mesh.elements.coordinate_transformation.QUAD_1d.X(quad_degree, quad_type)
        Q3_Y = mesh.elements.coordinate_transformation.QUAD_1d.Y(quad_degree, quad_type)
        Q3_JM = mesh.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, quad_type)
        Q3_J00 = mesh.elements.coordinate_transformation.QUAD_1d.J00(quad_degree, quad_type)
        Q3_J01 = mesh.elements.coordinate_transformation.QUAD_1d.J01(quad_degree, quad_type)
        Q3_J10 = mesh.elements.coordinate_transformation.QUAD_1d.J10(quad_degree, quad_type)
        Q3_J11 = mesh.elements.coordinate_transformation.QUAD_1d.J11(quad_degree, quad_type)
        Q3_J = mesh.elements.coordinate_transformation.QUAD_1d.Jacobian(quad_degree, quad_type)
        Q3_M = mesh.elements.coordinate_transformation.QUAD_1d.metric(quad_degree, quad_type)
        Q3_MM = mesh.elements.coordinate_transformation.QUAD_1d.metric_matrix(quad_degree, quad_type)
        Q3_iJM = mesh.elements.coordinate_transformation.QUAD_1d.inverse_Jacobian_matrix(quad_degree, quad_type)
        Q3_iJ = mesh.elements.coordinate_transformation.QUAD_1d.inverse_Jacobian(quad_degree, quad_type)
        Q3_iMM = mesh.elements.coordinate_transformation.QUAD_1d.inverse_metric_matrix(quad_degree, quad_type)

        for i in mesh.elements:
            np.testing.assert_array_almost_equal(_mapping[i], Q3_mapping[i])
            np.testing.assert_array_almost_equal(_X[i], Q3_X[i])
            np.testing.assert_array_almost_equal(_Y[i], Q3_Y[i])
            for j in range(2):
                for k in range(2):
                    np.testing.assert_array_almost_equal(_JM[i][j][k], Q3_JM[i][j][k])
            np.testing.assert_array_almost_equal(_J00[i], Q3_J00[i])
            np.testing.assert_array_almost_equal(_J01[i], Q3_J01[i])
            np.testing.assert_array_almost_equal(_J10[i], Q3_J10[i])
            np.testing.assert_array_almost_equal(_J11[i], Q3_J11[i])
            np.testing.assert_array_almost_equal(_J[i], Q3_J[i])
            np.testing.assert_array_almost_equal(_M[i], Q3_M[i])
            np.testing.assert_array_almost_equal(_MM[i], Q3_MM[i])
            np.testing.assert_array_almost_equal(_iJM[i], Q3_iJM[i])
            np.testing.assert_array_almost_equal(_iJ[i], Q3_iJ[i])
            np.testing.assert_array_almost_equal(_iMM[i], Q3_iMM[i])

    return 1


def test_Mesh_NO4_mesh_trace_topology():
    """Unittests for the mesh."""
    if rAnk == mAster_rank:
        print("+++ [test_Mesh_NO4_mesh_trace_topology] ...... ", flush=True)

    MID = list(DomainInputFinder.___defined_DI___().keys())
    if rAnk == mAster_rank:
        __ = random.sample(range(0,len(MID)), 4)
        meshes = [MID[i] for i in __]
        II = random.randint(3,4) # [II, JJ] element layout
        JJ = random.randint(2,5) # [II, JJ] element layout
    else:
        meshes = None
        II, JJ = None, None
    II, JJ = cOmm.bcast([II, JJ], root=mAster_rank)
    meshes = cOmm.bcast(meshes, root=mAster_rank)

    for mid in meshes:
        # ... generate meshes ...
        if mid in ('crazy', 'crazy_periodic'):
            if rAnk == mAster_rank:
                c = random.uniform(0, 0.3)
            else:
                c = None
            c = cOmm.bcast(c, root=mAster_rank)
            mesh = MeshGenerator(mid, c=c)([II, JJ], EDM='debug')
        else:
            mesh = MeshGenerator(mid)([II, JJ], EDM='debug')

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

    mesh = MeshGenerator('cic')([3, 2], EDM='debug')
    MAP = mesh.trace.elements.map
    if 1 in MAP:
        assert MAP[1] == [1, 4, 5, 6]
        e = mesh.trace.elements[1]
        assert e.positions ==('0D', '1U')
        assert e.CHARACTERISTIC_position in e.positions
        assert str(e.CHARACTERISTIC_element) + e.CHARACTERISTIC_edge == e.CHARACTERISTIC_position
        assert e.IS.on_periodic_boundary is False
        assert e.IS.on_mesh_boundary is False
    if 17 in MAP:
        assert MAP[17] == [43, 45, 40, 46]
        e = mesh.trace.elements[45]
        assert e.positions ==('17D', '21U')
        assert e.CHARACTERISTIC_position in e.positions
        assert str(e.CHARACTERISTIC_element) + e.CHARACTERISTIC_edge == e.CHARACTERISTIC_position
        assert e.IS.on_periodic_boundary is False
        assert e.IS.on_mesh_boundary is False
        e = mesh.trace.elements[46]
        assert e.positions ==('17R', 'Down')
        assert e.CHARACTERISTIC_position == '17R'
        assert e.CHARACTERISTIC_position in e.positions
        assert str(e.CHARACTERISTIC_element) + e.CHARACTERISTIC_edge == e.CHARACTERISTIC_position
        assert e.IS.on_periodic_boundary is False
        assert e.IS.on_mesh_boundary
    if 33 in MAP:
        assert MAP[33] == [81, 82, 76, 83]
        e = mesh.trace.elements[81]
        assert e.positions ==('33U', 'Upper')
        assert e.CHARACTERISTIC_position == '33U'
        assert e.CHARACTERISTIC_position in e.positions
        assert str(e.CHARACTERISTIC_element) + e.CHARACTERISTIC_edge == e.CHARACTERISTIC_position
        assert e.IS.on_periodic_boundary is False
        assert e.IS.on_mesh_boundary

    if 27 in MAP: assert MAP[27] == [67, 68, 62, 69]
    if 28 in MAP: assert MAP[28] == [68, 70, 64, 71]
    if 29 in MAP: assert MAP[29] == [70, 72, 66, 73]

    return 1





if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_2d\__tests__\unittests\mesh.py
    test_Mesh_NO1_mesh_topology()
    # test_Mesh_NO2_mesh_coordinate_transformation()
    # test_Mesh_NO3_mesh_coordinate_transformation_QUAD()
    # test_Mesh_NO4_mesh_trace_topology()