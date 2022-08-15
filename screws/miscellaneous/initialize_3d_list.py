# -*- coding: utf-8 -*-

def initialize_3d_list(a, b, c):
    """

    :param a:
    :param b:
    :param c:
    :return:
    """
    lst = [[[None for _ in range(c)] for _ in range(b)] for _ in range(a)]
    return lst