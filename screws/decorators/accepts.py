# -*- coding: utf-8 -*-
"""Authors are gratefully appreciated.

decorators return new functions, therefore will change .__code__ and so on.

@author: Yi Zhang (collecting). Created on Thu Feb 22 22:37:03 2018
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, Netherlands

"""
import numpy as np
import functools
lru_cache = functools.lru_cache
from screws.decorators.timeit.timeit_1 import timeit1
from screws.decorators.timeit.timeit_2 import timeit2

from screws.decorators.memoize.memoize_1 import memoize1
from screws.decorators.memoize.memoize_2 import memoize2
from screws.decorators.memoize.memoize_3 import memoize3
from screws.decorators.memoize.memoize_4 import memoize4
from screws.decorators.memoize.memoize_5 import memoize5



def accepts(*types):
    """ Works for methods and functions.
    
    When the first argument is self (for method), use ``@accepts('self', ...)``.
    
    - ``@accepts(..., '', ...)``: '' means everything is OK
    - ``@accepts(..., None, ...)``: None means everything is OK
    - ``@accepts('self', 'int', 'float', ('int', 'ndarray'))``
    - ``@accepts('self', 'int', 'float', ['int', 'ndarray'])``
    
    ``('int', 'ndarray')`` or ``['int', 'ndarray']`` is the 'multiple acceptance' case. So the
    input can be one of these types.
    
    class names in str will be always ok.
    If use type classes: int, float, str, make sure they are before the str class name. For example:

    - ``(int, float, 'ndarray')`` is GOOD.
    - ``('ndarray', int, float)`` is BAD.
    
    .. warning::

        It does NOT work with ``@classmethod``.
    
    But it works with ``@staticmethod`` in the way:
    ::

        @staticmethod
        @accepts(('int', "float"), 'int')
        def staticMethodName(a, b):
            ...
    
    We have defined some acceptance keys:

     1. sparse_matrix: including all sparse matrix types in `scipy.sparse`,
             'csc_matrix', 'csr_matrix', 'coo_matrix', 'lil_matrix', 
             'dok_matrix', 'dia_matrix', 'bsr_matrix'.
     2. natural_number: 0, 1, 2, 3, ......
     3. positive_int: 1, 2, 3, 4, ......
     4. negative_int: -1, -2, -3, -4, ......
    
    .. warning::

        Following self-defined types can not be used in 'multiple acceptance' case:
        ``natural_number``;  ``positive_int``; ``negative_int``.
        
    .. seealso::

        variable ``_not_in_multiple_acceptance_``.

    when we have shape requirement for ith input, then we need have to add a 
    shape keyword to it like:

    - ``shape=(i,j,k,...)``: The ith input has to be of shape ``(i,j,k,...)``.
    - ``shape=(4,x)``: The input can only be of 2-dimensions, and the first dimension must be of depth 4.
    - ``shape=(4,x,...)``: The input can only be of 2 or more dimensions, and the first dimension must be of depth 4.
    - ``shape=(x,3,...)``: the input can only be of 2 or more dimensions, and the second dimension must be of depth 3.
        
    .. warning::

        Currently, '...' can only be put at last, so ``shape:(...,1,2)`` is not allowed.
    
    For check dimension requirement, we use ``ndim=x``.
    """    
    def check_accepts(f):
        def new_f(*args, **kwds):
            for i, (a, t) in enumerate(zip(args, types)):
                # we handle iterable accepts in tuple
                if isinstance(t, list): t = tuple(t)
                #______________________________________________________________________
                # we have a self-defined type: sparse_matrix, which basically
                # includes all sparse type in scipy.sparse.
                if t == 'sparse_matrix':
                    assert a.__class__.__name__ in ('csc_matrix', 'csr_matrix', 'coo_matrix', 'lil_matrix', 
                                                    'dok_matrix', 'dia_matrix', 'bsr_matrix'), \
                        " <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t)
                elif t == 'natural_number':
                    assert isinstance(a, int) and a >= 0, \
                        " <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t)
                elif t == 'positive_int':
                    assert isinstance(a, int) and a > 0, \
                        " <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t)
                elif t == 'negative_int':
                    assert isinstance(a, int) and a < 0, \
                        " <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t)
                elif isinstance(t, tuple):
                    if 'sparse_matrix' in t:
                        t += ('csc_matrix', 'csr_matrix', 'coo_matrix', 'lil_matrix', 
                              'dok_matrix', 'dia_matrix', 'bsr_matrix')
                    else:
                        pass
                    _not_in_multiple_acceptance_ = ('natural_number', 'positive_int', 'negative_int')
                    assert all([ti not in _not_in_multiple_acceptance_ for ti in t]), \
                    " <{}> : {}th Arg: {} has not-multiple-acceptance type.".format(f.__name__, i, t)
                else:
                    pass
                #______________________________________________________________________
                if t is None or t == '' or (t == 'self' and i == 0):
                    # When used for method, we set the first to be 'self' of course.
                    # When t == '' or is None, then the input can be anything.
                    pass
                elif t in ('sparse_matrix', 'natural_number', 'positive_int', 'negative_int'):
                    # those types already are judged.
                    pass
                else:
                    # shape requirement_________________________________________________
                    if isinstance(t, (list, tuple)):
                        for i_type in t:
                            #~~~~~~~~~~ check shape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            if isinstance(i_type, str) and i_type[0:5] == 'shape':
                                input_shape = np.shape(a)
                                input_ndim = np.ndim(a)
                                if input_shape == ():
                                    raise Exception(
                                    " <{}> : {}th Arg has not shape, so can not match shape requirement.".format(
                                    f.__name__, i))
                                assert '=' in i_type
                                shape = i_type.split('=')[1]
                                assert shape[0] == '(' and shape[-1] == ')', \
                                    " <{}> : {}th Arg shape requirement={} is wrong.".format(
                                                        f.__name__, i, i_type)
                                shape = shape[1:-1].split(',')
                                ls = len(shape)
                                #_____Future potential changes below___________________
                                for j, sj in enumerate(shape):
                                    if sj == '...':
                                        assert j+1 == ls, \
                                        " <{}> : for shape requirement, '...' must be put at last".format(
                                                f.__name__)
                                    elif sj == 'x':
                                        pass
                                    else:
                                        try:
                                            # noinspection PyTypeChecker
                                            shape[j] = int(sj)
                                        except ValueError:
                                            raise Exception(
                                                " <{}> : {}th Arg shape requirement={} is wrong.".format(
                                                        f.__name__, i, i_type))
                                require_shape = tuple(shape)
                                try:
                                    for j, rsk in enumerate(require_shape):
                                        if rsk == '...':
                                            break
                                        elif rsk == 'x':
                                            pass
                                        else:
                                            assert input_ndim > j and input_shape[j] == rsk
                                except AssertionError:
                                    raise Exception(
                                        " <{}> : {}th Arg does not match shape requirement: {}.".format(
                                        f.__name__, i, require_shape))
                                #------------------------------------------------------
                                break
                            #~~~~~~~~~~ check ndim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            elif isinstance(i_type, str) and i_type[0:4] == 'ndim':
                                assert '=' in i_type
                                ndim = int(i_type.split('=')[1])
                                assert np.ndim(a) == ndim, \
                                " <{}> : {}th Arg does not match ndim requirement: {}.".format(
                                        f.__name__, i, ndim)
                            #~~~~~~~~~~ ELSE: pass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            else:
                                pass
                            #==========================================================
                    else:
                        pass
                    # type requirement__________________________________________________
                    try:
                        assert isinstance(a, t)
                    except TypeError:
                        try:
                            if isinstance(t, str):
                                t = (t,)
                            assert a.__class__.__name__ in t
                        except TypeError:
                            assert a.__class__.__name__ == t, " <{}> : {}th Arg: %r does not match: %s".format(
                                    f.__name__, i)%(a,t)
                        except AssertionError:
                            raise Exception(" <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t))
                    except AssertionError:
                        raise Exception(" <{}> : {}th Arg: %r does not match: %s".format(f.__name__, i)%(a,t))
                    #------------------------------------------------------------------
                #----------------------------------------------------------------------
            return f(*args, **kwds)
        new_f.__name__ = f.__name__
        return new_f
    return check_accepts



















if __name__ == "__main__":
    # do some tests
    @timeit2
    def network_call(user_id):
        print("(computed)")
        return user_id

    class NetworkEngine(object):
        def __init__(self):
            pass

        @lru_cache()  # test for #1,2
        def search(self, user_id):
            return network_call(user_id)

        @staticmethod
        def test_tt(user_id):
            return network_call(user_id)

    e = NetworkEngine()
    for v in [1,2,3,3,3,1,4,1,5,'abs','abs','abs']:
        print(e.search(v))

    @accepts('positive_int', (list, float, list, 'ndim=1', 'shape=(3)'))
    def foo(a, b):
        return a, b
    
    foo(1, [1,2,3])

    @memoize4
    @timeit1
    def foo1(a, b, c=-2):
        print('computed')
        return a * b * c

    print(foo1(3, 5))
    print(foo1(4, 5))
    print(foo1(5, 5))
    print(foo1(5, 5))
    print(foo1(4, 5))
    
    ms = memoize1, memoize2, memoize3, memoize4

    @memoize1
    def fun2(a):
        print(a)
        

    @memoize5
    def fun5(b):
        print(b)