

import numpy as np

class RegionTopology:
    @classmethod
    def num_sides(cls):
        return 6

    @classmethod
    def num_edges(cls):
        return 12

    @classmethod
    def num_corners(cls):
        return 8

    @classmethod
    def _side_name_to_index_(cls, _side_name_):
        return {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}[_side_name_]

    @classmethod
    def _side_index_to_name_(cls, _side_index_):
        return {0: 'N', 1: 'S', 2: 'W', 3: 'E', 4: 'B', 5: 'F'}[_side_index_]

    @classmethod
    def _side_name_to_index_dict_(cls):
        return {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}

    @classmethod
    def _side_index_to_name_dict_(cls):
        return {0: 'N', 1: 'S', 2: 'W', 3: 'E', 4: 'B', 5: 'F'}

    @classmethod
    def _side_pairing_(cls):
        return {'N': 'S', 'S': 'N', 'W': 'E', 'E': 'W', 'B': 'F', 'F': 'B'}

    @classmethod
    def _side_index_pairing_(cls):
        return {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}

    @classmethod
    def _corner_name_to_index_(cls, _corner_name_):
        return {'NWB': 0, 'SWB': 1, 'NEB': 2, 'SEB': 3, 'NWF': 4, 'SWF': 5, "NEF": 6, 'SEF': 7,
                'NBW': 0, 'SBW': 1, 'NBE': 2, 'SBE': 3, 'NFW': 4, 'SFW': 5, "NFE": 6, 'SFE': 7,
                'WNB': 0, 'WSB': 1, 'ENB': 2, 'ESB': 3, 'WNF': 4, 'WSF': 5, "ENF": 6, 'ESF': 7,
                'WBN': 0, 'WBS': 1, 'EBN': 2, 'EBS': 3, 'WFN': 4, 'WFS': 5, "EFN": 6, 'EFS': 7,
                'BNW': 0, 'BSW': 1, 'BNE': 2, 'BSE': 3, 'FNW': 4, 'FSW': 5, "FNE": 6, 'FSE': 7,
                'BWN': 0, 'BWS': 1, 'BEN': 2, 'BES': 3, 'FWN': 4, 'FWS': 5, "FEN": 6, 'FES': 7}[_corner_name_]

    @classmethod
    def _corner_index_to_name_(cls, _corner_index_):
        return {0: 'NWB', 1: 'SWB', 2: 'NEB', 3: 'SEB', 4: 'NWF', 5: 'SWF', 6: "NEF", 7: 'SEF'}[_corner_index_]

    @classmethod
    def _corner_name_to_index_dict_(cls):
        return {'NWB': 0, 'SWB': 1, 'NEB': 2, 'SEB': 3, 'NWF': 4, 'SWF': 5, "NEF": 6, 'SEF': 7,
                'NBW': 0, 'SBW': 1, 'NBE': 2, 'SBE': 3, 'NFW': 4, 'SFW': 5, "NFE": 6, 'SFE': 7,
                'WNB': 0, 'WSB': 1, 'ENB': 2, 'ESB': 3, 'WNF': 4, 'WSF': 5, "ENF": 6, 'ESF': 7,
                'WBN': 0, 'WBS': 1, 'EBN': 2, 'EBS': 3, 'WFN': 4, 'WFS': 5, "EFN": 6, 'EFS': 7,
                'BNW': 0, 'BSW': 1, 'BNE': 2, 'BSE': 3, 'FNW': 4, 'FSW': 5, "FNE": 6, 'FSE': 7,
                'BWN': 0, 'BWS': 1, 'BEN': 2, 'BES': 3, 'FWN': 4, 'FWS': 5, "FEN": 6, 'FES': 7}

    @classmethod
    def _corner_index_to_name_dict_(cls):
        return {0: 'NWB', 1: 'SWB', 2: 'NEB', 3: 'SEB', 4: 'NWF', 5: 'SWF', 6: "NEF", 7: 'SEF'}

    @classmethod
    def _side_corner_local_numbering_(cls, _side_index_):
        """
        Notice the values are always increasing. The reason is the same with the reason
        for `_side_edge_local_numbering_` or `_side_edge_local_numbering_dict_`.

        """
        return {0: (0, 2, 4, 6),  # the north side has corners locally numbered 0, 2, 4, 6.
                1: (1, 3, 5, 7),  # the south side has corners locally numbered 1, 3, 5, 7.
                2: (0, 1, 4, 5),  # the west side has corners locally numbered 0, 1, 4 5.
                3: (2, 3, 6, 7),  # and so on
                4: (0, 1, 2, 3),
                5: (4, 5, 6, 7)}[_side_index_]

    @classmethod
    def _side_corner_local_numbering_dict_(cls):
        """
        Notice the values are always increasing. The reason is the same with the reason
        for `_side_edge_local_numbering_` or `_side_edge_local_numbering_dict_`.

        """
        return {0: (0, 2, 4, 6),  # the north side has corners locally numbered 0, 2, 4, 6.
                1: (1, 3, 5, 7),  # the south side has corners locally numbered 1, 3, 5, 7.
                2: (0, 1, 4, 5),  # the west side has corners locally numbered 0, 1, 4 5.
                3: (2, 3, 6, 7),  # and so on
                4: (0, 1, 2, 3),
                5: (4, 5, 6, 7)}

    @classmethod
    def _side_edge_local_numbering_(cls, _side_index_):
        """
        Notice the values are always increasing. This is because we are consistently
        using the same convention of numbering local geometries. dx goes first, dy then
        , dz the last.

        """
        return {0: (4, 6, 8, 10),  # the north side has edges locally numbered 4,6,8,10
                1: (5, 7, 9, 11),  # the south side has edges locally numbered 5,7,9,11
                2: (0, 2, 8, 9),  # the west side has edges locally numbered 0,2,8,9
                3: (1, 3, 10, 11),  # and so on
                4: (0, 1, 4, 5),
                5: (2, 3, 6, 7)}[_side_index_]

    @classmethod
    def _side_edge_local_numbering_dict_(cls):
        """
        Notice the values are always increasing. This is because we are consistently
        using the same convention of numbering local geometries. dx goes first, dy then
        , dz the last.

        """
        return {0: (4, 6, 8, 10),  # the north side has edges locally numbered 4,6,8,10
                1: (5, 7, 9, 11),  # the south side has edges locally numbered 5,7,9,11
                2: (0, 2, 8, 9),  # the west side has edges locally numbered 0,2,8,9
                3: (1, 3, 10, 11),  # and so on
                4: (0, 1, 4, 5),
                5: (2, 3, 6, 7)}

    @classmethod
    def _side_axis_distribution_(cls):
        """
        Here 'N':(0,0) means the N side is perpendicular to the 0-axis and at
        the starting side. 'S':(0,-1) means the 'D' side is perpendicular to
        the 0-axis but at the end side, so -1. And 'F':(2,-1) means the F side
        is perpendicular to 2-axis and at the end side.

        And so on.

        """
        return {'N': (0, 0), 'S': (0, -1), 'W': (1, 0), 'E': (1, -1), 'B': (2, 0), 'F': (2, -1)}

    @classmethod
    def _axis_indix_dict_(cls):
        return {'x': 0, 'y': 1, 'z': 2}

    @classmethod
    def _edge_name_to_index_(cls, _name_):
        return \
            {'WB': 0, 'EB': 1, 'WF': 2, 'EF': 3, 'NB': 4, 'SB': 5, 'NF': 6, 'SF': 7, 'NW': 8, 'SW': 9, 'NE': 10,
             'SE': 11,
             'BW': 0, 'BE': 1, 'FW': 2, 'FE': 3, 'BN': 4, 'BS': 5, 'FN': 6, 'FS': 7, 'WN': 8, 'WS': 9, 'EN': 10,
             'ES': 11}[
                _name_]

    @classmethod
    def _edge_index_to_name_(cls, _index_):
        return \
            {0: 'WB', 1: 'EB', 2: 'WF', 3: 'EF', 4: 'NB', 5: 'SB', 6: 'NF', 7: 'SF', 8: 'NW', 9: 'SW', 10: 'NE',
             11: 'SE'}[
                _index_]

    @classmethod
    def _edge_name_to_index_dict_(cls):
        return {'WB': 0, 'EB': 1, 'WF': 2, 'EF': 3, 'NB': 4, 'SB': 5, 'NF': 6, 'SF': 7, 'NW': 8, 'SW': 9, 'NE': 10,
                'SE': 11,
                'BW': 0, 'BE': 1, 'FW': 2, 'FE': 3, 'BN': 4, 'BS': 5, 'FN': 6, 'FS': 7, 'WN': 8, 'WS': 9, 'EN': 10,
                'ES': 11}

    @classmethod
    def _edge_index_to_name_dict_(cls):
        return {0: 'WB', 1: 'EB', 2: 'WF', 3: 'EF', 4: 'NB', 5: 'SB', 6: 'NF', 7: 'SF', 8: 'NW', 9: 'SW', 10: 'NE',
                11: 'SE'}

    @classmethod
    def _corner_local_numbering_(cls):
        """ The local numbering of corners of 3D regions. """
        return np.array(range(8)).reshape((2, 2, 2), order='F')

    @classmethod
    def _edge_local_numbering_(cls):
        return (np.arange(4).reshape((1, 2, 2), order='F'),
                np.arange(4).reshape((2, 1, 2), order='F') + 4,
                np.arange(4).reshape((2, 2, 1), order='F') + 8)

    @classmethod
    def _side_local_numbering_(cls):
        return (np.arange(2).reshape((2, 1, 1), order='F'),
                np.arange(2).reshape((1, 2, 1), order='F') + 2,
                np.arange(2).reshape((1, 1, 2), order='F') + 4)