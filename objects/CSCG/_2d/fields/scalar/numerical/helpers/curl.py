# -*- coding: utf-8 -*-
"""

"""


class ___VECTOR_CURL_HELPER___(object):
    def __init__(self, f):
        self._f_ = f

    def __call__(self, t, x, y):
        return - self._f_(t, x, y)
