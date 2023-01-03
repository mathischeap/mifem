# -*- coding: utf-8 -*-
from components.decorators.all import accepts


@accepts((int, float), (int, float), (int, float))
def check_almost_in_range(a, lb, ub, tol=1e-8):
    """check a is in [lb-tol, up+tol].

    :param a:
    :param lb:
    :param ub:
    :param tol:
    :return:
    """
    assert lb < ub, f"lower bound {lb} must be lower than upper bound {ub}."
    return (lb - tol) < a < (ub + tol)
