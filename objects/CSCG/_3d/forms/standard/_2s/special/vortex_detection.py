# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
import numpy as np


class ___3dCSCG_2Form_Vortex_Detection___(FrozenOnly):
    """A wrapper of all vortex detection methods. So, we consider this 1 form as
    a variable of a flow field."""
    def __init__(self, _2sf):
        self._sf_ = _2sf
        self._freeze_self_()

    def ___PRIVATE_generate_gradient_tensor_at___(self, xi, eta, sigma):
        """We compute the gradient tensor of this 2form. To do so, we first project
        this 2-form into a vector of 3 standard 0-forms which represent the three
        components. Then we do the gradient (apply the incidence matrix E10) to each
        standard 0-form.

        It returns a 3 by 3 tensor representing
            ((du_dx, du_dy, du_dz),
             (dv_dx, dv_dy, dv_dz),
             (dw_dx, dw_dy, dw_dz)).
        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]

        """
        assert np.ndim(xi) == 1 and np.all(np.diff(xi) >0) and np.max(xi) <= 1 and np.min(xi) >= -1, \
            f"xi={xi} wrong, should be 1d array in [-1,1] and increasing."
        assert np.ndim(eta) == 1 and np.all(np.diff(eta) >0) and np.max(eta) <= 1 and np.min(eta) >= -1, \
            f"eta={eta} wrong, should be 1d array in [-1,1] and increasing."
        assert np.ndim(sigma) == 1 and np.all(np.diff(sigma) >0) and np.max(sigma) <= 1 and np.min(sigma) >= -1, \
            f"sigma={sigma} wrong, should be 1d array in [-1,1] and increasing."

        U, V, W = self._sf_.projection.to.vector_of_3_standard_0forms()
        dU = U.coboundary()
        dV = V.coboundary()
        dW = W.coboundary()

        xyz, rdU = dU.reconstruct(xi, eta, sigma, ravel=False)
        xyz, rdV = dV.reconstruct(xi, eta, sigma, ravel=False)
        xyz, rdW = dW.reconstruct(xi, eta, sigma, ravel=False)

        dU_dx, dU_dy, dU_dz = dict(), dict(), dict()
        dV_dx, dV_dy, dV_dz = dict(), dict(), dict()
        dW_dx, dW_dy, dW_dz = dict(), dict(), dict()

        for i in rdU:
            dU_dx[i], dU_dy[i], dU_dz[i] = rdU[i]
            dV_dx[i], dV_dy[i], dV_dz[i] = rdV[i]
            dW_dx[i], dW_dy[i], dW_dz[i] = rdW[i]

        return xyz, ((dU_dx, dU_dy, dU_dz),
                     (dV_dx, dV_dy, dV_dz),
                     (dW_dx, dW_dy, dW_dz))

    def ___PRIVATE_generate_S_and_Omega___(self, xi, eta, sigma):
        """
        S and Omega are the symmetric and antisymmetric components of gradient
        tensor G. So both S and Omega are 3 by 3 tensor, and

        S_{i,j} = 0.5 * (G_{i,j} + G_{j,i})
        Omega_{i,j} = 0.5 * (G_{i,j} - G_{j,i})

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]

        """
        S_00, S_01, S_02 = dict(), dict(), dict()
        S_10, S_11, S_12 = dict(), dict(), dict()
        S_20, S_21, S_22 = dict(), dict(), dict()

        O_00, O_01, O_02 = dict(), dict(), dict()
        O_10, O_11, O_12 = dict(), dict(), dict()
        O_20, O_21, O_22 = dict(), dict(), dict()

        xyz, GT = self.___PRIVATE_generate_gradient_tensor_at___(xi, eta, sigma)

        dU_xyz, dV_xyz, dW_xyz = GT
        U_00, U_01, U_02 = dU_xyz
        U_10, U_11, U_12 = dV_xyz
        U_20, U_21, U_22 = dW_xyz

        for i in U_00: # will go through all local mesh elements
            u_00 = U_00[i]
            u_01 = U_01[i]
            u_02 = U_02[i]
            u_10 = U_10[i]
            u_11 = U_11[i]
            u_12 = U_12[i]
            u_20 = U_20[i]
            u_21 = U_21[i]
            u_22 = U_22[i]

            s_00 = u_00
            s_11 = u_11
            s_22 = u_22

            s_01 = 0.5 * (u_01 + u_10)
            s_02 = 0.5 * (u_02 + u_20)
            s_10 = s_01
            s_12 = 0.5 * (u_12 + u_21)
            s_20 = s_02
            s_21 = s_12
            S_00[i], S_01[i], S_02[i] = s_00, s_01, s_02
            S_10[i], S_11[i], S_12[i] = s_10, s_11, s_12
            S_20[i], S_21[i], S_22[i] = s_20, s_21, s_22

            o_00 = o_11 = o_22 = 0
            o_01 = 0.5 * (u_01 - u_10)
            o_02 = 0.5 * (u_02 - u_20)
            o_10 = - o_01
            o_12 = 0.5 * (u_12 - u_21)
            o_20 = - o_02
            o_21 = - o_12
            O_00[i], O_01[i], O_02[i] = o_00, o_01, o_02
            O_10[i], O_11[i], O_12[i] = o_10, o_11, o_12
            O_20[i], O_21[i], O_22[i] = o_20, o_21, o_22

        return xyz, ((S_00, S_01, S_02), (S_10, S_11, S_12), (S_20, S_21, S_22)), \
                    ((O_00, O_01, O_02), (O_10, O_11, O_12), (O_20, O_21, O_22))

    def ___PRIVATE_generate_lambda_1_2_3_Q___(self, xi, eta, sigma):
        """  See [on the identification of a vortex] by Jeong and Hussain.

        Compute the lambda_2, and Q definitions.

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]
        """
        xyz, S, O = self.___PRIVATE_generate_S_and_Omega___(xi, eta, sigma)

        S0, S1, S2 = S
        O0, O1, O2 = O

        S00, S01, S02 = S0
        S10, S11, S12 = S1
        S20, S21, S22 = S2

        O00, O01, O02 = O0
        O10, O11, O12 = O1
        O20, O21, O22 = O2

        Q = dict()
        LAMBDA_2 = dict()

        for i in S00: # we go through all local mesh elements
            s00, s01, s02 = S00[i], S01[i], S02[i]
            s10, s11, s12 = S10[i], S11[i], S12[i]
            s20, s21, s22 = S20[i], S21[i], S22[i]
            o00 ,o01, o02 = O00[i], O01[i], O02[i]
            o10, o11, o12 = O10[i], O11[i], O12[i]
            o20, o21, o22 = O20[i], O21[i], O22[i]

            SHAPE = s00.shape

            so_00 = (s00**2 + o00**2).ravel('F')
            so_01 = (s01**2 + o01**2).ravel('F')
            so_02 = (s02**2 + o02**2).ravel('F')

            so_10 = (s10**2 + o10**2).ravel('F')
            so_11 = (s11**2 + o11**2).ravel('F')
            so_12 = (s12**2 + o12**2).ravel('F')

            so_20 = (s20**2 + o20**2).ravel('F')
            so_21 = (s21**2 + o21**2).ravel('F')
            so_22 = (s22**2 + o22**2).ravel('F')

            so = np.zeros((len(so_00),3,3))
            so[:,0,0] = so_00
            so[:,0,1] = so_01
            so[:,0,2] = so_02
            so[:,1,0] = so_10
            so[:,1,1] = so_11
            so[:,1,2] = so_12
            so[:,2,0] = so_20
            so[:,2,1] = so_21
            so[:,2,2] = so_22

            eigen_values, _ = np.linalg.eig(so)
            Q[i] = (- 0.5 * np.sum(eigen_values, axis=1)).reshape(SHAPE, order='F')
            eigen_values = np.sort(eigen_values, axis=1)

            # print(eigen_values, eigen_values[:,1] )

            LAMBDA_2[i] = eigen_values[:,1].reshape(SHAPE, order='F')

        return xyz, Q, LAMBDA_2

    def Q_and_lambda2(self, xi, eta, sigma):
        """ See [on the identification of a vortex] by Jeong and Hussain.

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]
        """
        xyz, Q, LAMBDA_2 = self.___PRIVATE_generate_lambda_1_2_3_Q___(xi, eta, sigma)[:3]
        return xyz, Q, LAMBDA_2

