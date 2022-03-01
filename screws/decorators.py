# -*- coding: utf-8 -*-
"""
Collections of some memorizers. These memorizers are found on internet.
Authors are gratefully appreciated.

decorators return new functions, therefore will change .__code__ and so on.

@author: Yi Zhang (collecting). Created on Thu Feb 22 22:37:03 2018
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, Netherlands
         
"""
import numpy as np
import functools
from time import localtime, strftime, time
lru_cache = functools.lru_cache


def memoize1(func):
    """
    Generally speaking, the best one is this one.

    - ``+``: Can be used for frozen object.
    - ``+``: Can be used for numpy.ndarray inputs. But normally when we can have numpy.ndarray as input,
        we do not use @memoize because storing the input may need a lot of memory.
    - ``-``: it is relatively slower than others.
    - ``-``: kwargs seem not to cached in keys at all.
    """
    cache = func.cache = {}

    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)
        try:
            return cache[key]
        except KeyError:
            cache[key] = func(*args, **kwargs)
            return cache[key]

    return memoized_func


class memoize2(object):
    """
    Cache the return value of a method in class.

    This class is meant to be used as a decorator of methods. The return value from
    a given method invocation will be cached on the instance whose method was
    invoked. All arguments passed to a method decorated with memoize must be
    hashable.

    If a memoized method is invoked directly on its class the result will not be
    cached. Instead the method will be invoked like a static method:
    ::

       class Obj(object):
           @memoize2
           def add_to(self, arg):
               return self + arg

    ``Obj.add_to(1)`` # not enough arguments; ``Obj.add_to(1, 2)`` # returns 3, result is not cached

    Script borrowed from here:
    MIT Licensed, attributed to Daniel Miller, Wed, 3 Nov 2010
    `active-state <http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/>`_

    - ``-``: Can not be used for frozen object. To use it for frozen object, use the method once before
        freeze the object.
    - ``-``: Can not be used for numpy.ndarray inputs.
    - ``+``: faster than memoize1.

    - ``-``: can not be dumped when used together with @accepts for one method.
    """

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return functools.partial(self, obj)

    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res


def memoize3(f):
    """
    Memoization decorator for a function taking a single argument.

    - ``+``: Very fast.
    - ``-``: Can not be used for methods in classes.
    - ``-``: for single input functions.
    - ``-``: Can not be used for numpy.ndarray inputs.
    """

    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret

    return memodict().__getitem__


def memoize4(f):
    """
    Memoization decorator for a function taking one or more arguments.

    - ``+``: Very fast.
    - ``+``: Can be used for multiple inputs functions.
    - ``-``: Can not be used for methods in classes.
    - ``-``: Can not be used for numpy.ndarray inputs.
    """

    class memodict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret

    return memodict().__getitem__


class memoize5(object):
    """
    cache the return value of a method.

    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.

    If a memoized method is invoked directly on its class the result will not
    be cached. Instead of the method will be invoked like a static method:
    class Obj(object):
        @memoize
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    """

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return functools.partial(self, obj)

    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res

# from wrapt import synchronized
# synchronized
# + : this decorator works for functions, classes, and methods. 
# - : It will create new attributes when used for methods, so make sure the class is melt or use it at least once
#     before freeze the class or manually define an attribute called _synchronized_lock before.
# - : when use this decorator, the object can not be pickled, so can not be saved.
# synchronized decorator is used to prohibit multiple threads are executing
# the same one method at the same time. This method usually produces attributes
# for the class. So if it is called at the same time, there may be confusions.


# from numba import jit

def accepts(*types):
    """ 
    works for methods and functions.
    
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



def timeit1(method):
    """A timer decorator for functions or methods."""

    def timed(*args, **kwargs):
        print(" <TimeIt> : Method [%r] with inputs: (%r, %r)" % (method.__name__, args, kwargs))
        print("            started at " + strftime("%Y-%m-%d %H:%M:%S", localtime()))
        ts = time()
        result = method(*args, **kwargs)
        minutes, seconds = divmod(time() - ts, 60)
        hours, minutes = divmod(minutes, 60)
        print(" <TimeIt> : Method [%r] Done, costs： %d:%02d:%02d (hh:mm:ss)" % (
            method.__name__, hours, minutes, seconds))
        return result

    return timed


def timeit2(method):
    """A timer decorator for functions or methods."""

    def timed(*args, **kwargs):
        print(" <TimeIt> -- Method [%r]:" %method.__name__)
        print("          -> Starts at " + strftime("%Y-%m-%d %H:%M:%S", localtime()))
        ts = time()
        result = method(*args, **kwargs)
        minutes, seconds = divmod(time() - ts, 60)
        hours, minutes = divmod(minutes, 60)
        print("          -> Done, costs %d:%02d:%02d (hh:mm:ss)" % (hours, minutes, seconds))
        return result

    return timed



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



if __name__ == "__main__":
    e = NetworkEngine()
    for v in [1,2,3,3,3,1,4,1,5,'abs','abs','abs']:
        print(e.search(v))

#    from scipy import sparse
#    a = sparse.lil_matrix((5, 5))

    @accepts('positive_int', (list, float, list, 'ndim=1', 'shape=(3)'))
    def foo(a, b):
        return a, b
    
    foo(1, [1,2,3])

    @memoize4
    @timeit1
#    @lru_cache(maxsize=2)
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
