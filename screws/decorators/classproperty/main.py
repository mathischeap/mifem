# -*- coding: utf-8 -*-

from screws.decorators.classproperty.descriptor import ClassPropertyDescriptor



def classproperty(func):
    if not isinstance(func, (classmethod, staticmethod)):
        func = classmethod(func)

    return ClassPropertyDescriptor(func)
