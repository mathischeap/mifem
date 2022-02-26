



import warnings


class MyCustomizedWarning(UserWarning, ValueError):
    pass


warnings.warn("deprecated", DeprecationWarning)



warnings.warn("My Customized Warning", MyCustomizedWarning)