"""

"""



from SCREWS.frozen import FrozenOnly
from importlib import import_module




class PreconditionerAllocator(FrozenOnly):
    """"""
    def __init__(self, ID):
        assert ID in self.___defined_preconditioners___(), f"Preconditioner ID={ID} is not implemented yet. Please use " \
                                                           f"one of {self.___defined_preconditioners___().keys()}"
        cls_name = self.___defined_preconditioners___()[ID]
        cls_path = self.___preconditioners_path___() + ID
        self._preconditioner_class_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, A, **kwargs):
        """"""
        assert str(A.__class__) == "<class 'TOOLS.linear_algebra.data_structures.GlobalMatrix'>", \
            f"A needs to be a 'TOOLS.linear_algebra.data_structures.GlobalMatrix'. Now I get {A.__class__}."
        return self._preconditioner_class_(A, **kwargs)

    @classmethod
    def ___defined_preconditioners___(cls):
        """Here we store all defined meshComponents. Whenever we define a new meshComponents (actually, a new
        domain_input), we add a nickname for it here.

        """
        _dict_ = {'Jacobian': "JacobianPreconditioner",
                 }
        return _dict_

    @classmethod
    def ___preconditioners_path___(cls):
        """ """
        return 'TOOLS.linear_algebra.preconditioners.'