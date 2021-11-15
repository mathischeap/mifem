"""
mapping a surface [xi, et] into (x,y,z)

xi = [-1, 1]
et = [-1, 1]


x = xi + 0.1*sin(2*pi*xi)sin(2*pi*et)
y = et + 0.1*sin(2*pi*xi)sin(2*pi*et)
z = 0.1*sin(2*pi*xi)sin(2*pi*et)

"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from numpy import sin, cos, pi

def X(xi, et):
    return xi + 0.1*sin(pi*xi) * sin(pi*et)

def Y(xi, et):
    return et + 0.1*sin(pi*xi) * sin(pi*et)

def Z(xi, et):
    return 0.1*sin(pi*xi) * sin(pi*et)


xi = np.linspace(-1,1,1000)
et = np.linspace(-1,1,1000)
xi, et = np.meshgrid(xi, et)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.view_init(45,60)
# ax.plot_surface(X(xi,et), Y(xi,et), Z(xi,et))
# plt.show()


def X_xi(xi, et):
    return 1 + 0.1*pi*cos(pi*xi) * sin(pi*et)

def X_et(xi, et):
    return 0.1*pi*sin(pi*xi) * cos(pi*et)


def Y_xi(xi, et):
    return 0.1*pi*cos(pi*xi) * sin(pi*et)

def Y_et(xi, et):
    return 1 + 0.1*pi*sin(pi*xi) * cos(pi*et)


def Z_xi(xi, et):
    return 0.1*pi*cos(pi*xi) * sin(pi*et)

def Z_et(xi, et):
    return 0.1*pi*sin(pi*xi) * cos(pi*et)


def metric_matrix(xi, et):
    J11 = X_xi(xi, et)
    J12 = X_et(xi, et)
    J21 = Y_xi(xi, et)
    J22 = Y_et(xi, et)
    J31 = Z_xi(xi, et)
    J32 = Z_et(xi, et)

    g_11 = J11**2 + J21**2 + J31**2
    g_22 = J12**2 + J22**2 + J32**2
    g_12 = J11*J12 + J21*J22 + J31*J32
    g_21 = g_12
    return ((g_11, g_12),
            (g_21, g_22))


def metric(xi, et):
    G = metric_matrix(xi, et)
    g_11 = G[0][0]
    g_12 = G[0][1]
    g_21 = g_12
    g_22 = G[1][1]
    return g_11 * g_22 - g_12 * g_21


g = metric(xi, et)


def area_element(xi, et):
    a1 = X_xi(xi, et)
    a2 = Y_xi(xi, et)
    a3 = Z_xi(xi, et)

    b1 = X_et(xi, et)
    b2 = Y_et(xi, et)
    b3 = Z_et(xi, et)

    c1 = a2*b3 - a3*b2
    c2 = a3*b1 - a1*b3
    c3 = a1*b2 - a2*b1

    return np.sqrt(c1**2+c2**2+c3**2)


AE = area_element(xi, et)
print(np.max(np.abs(np.sqrt(g)-AE)))

