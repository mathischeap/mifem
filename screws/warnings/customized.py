"""
"""





class MyCustomizedWarning(UserWarning, ValueError):
    pass



# import warnings
# warnings.warn("deprecated", DeprecationWarning)
#
# warnings.warn("My Customized Warning", MyCustomizedWarning)