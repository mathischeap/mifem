# -*- coding: utf-8 -*-
"""In this script, we DO NOT use the structure of naming files and folders of the mifem library.

"""



class StatisticError(Exception):
    """Raise when we try to access standard property ``statistic`` but ``___statistic___`` is not defined."""

class ParametersError(Exception):
    """Raise when we try to access standard property ``parameters`` but ``___parameters___`` is not defined."""

class FrozenError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class MeshError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class ElementsLayoutError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class ElementSidePairError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class ElementEdgePairError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class DimensionError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class LocalCochainShapeError(Exception):
    """Raise when we try to define new attribute for a frozen object."""

class InternetDisconnectedError(Exception):
    """Raise when we try to define new attribute for a frozen object."""


class LinerSystemSolverDivergenceError(Exception):
    """Raise when we try to define new attribute for a frozen object."""


class ThreeDimensionalTransfiniteInterpolationError(Exception):
    """Raise when we try to define new attribute for a frozen object."""


class EmailSendingError(Exception):
    """"""

class _2nCSCG_RF2_SegmentComparisonError(Exception):
    """"""