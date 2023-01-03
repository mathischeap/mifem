# -*- coding: utf-8 -*-


from components.decorators.all import accepts


@accepts((int, float), (int, float))
def check_multiple_close(a, b, tol=1e-8):
    """check if a = b*i +- tol where i = 1,2,3,4,...

    :param a:
    :param b:
    :param tol:
    :return:
    """
    remainder = a % b
    if remainder < tol:
        return True
    else:
        assert b > remainder, "something wrong."
        if (b - remainder) < tol:
            return True
    return False
