"""
Here we store some customize generators.

"""


from SCREWS.frozen import FrozenOnly



class CustomizedGenerator(FrozenOnly):
    """Not completed."""
    def __init__(self, data_generator, RANGE):
        """ """
        self._DG_ = data_generator
        self._range_ = RANGE
        self._freeze_self_()


    def __getitem__(self, i):
        return self._DG_(i)

    def __contains__(self, i):
        return i in self._range_

    def __iter__(self):
        for i in self._range_:
            yield i

    def __len__(self):
        return len(self._range_)