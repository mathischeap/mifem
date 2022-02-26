





class MyCustomizedWarning(UserWarning, ValueError):
    pass



class TraceElementWarning(UserWarning, ValueError):
    pass


# import warnings
# warnings.warn("deprecated", DeprecationWarning)
#
# warnings.warn("My Customized Warning", MyCustomizedWarning)