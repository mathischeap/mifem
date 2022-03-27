



class RegionTopology:
    @classmethod
    def num_edges(cls):
        return 4

    @classmethod
    def num_corners(cls):
        return 4

    @classmethod
    def _edge_name_to_index_dict_(cls):
        return {'U': 0, 'D': 1, 'L': 2, 'R': 3}

    @classmethod
    def _edge_name_to_index_(cls, _edge_name_):
        return {'U': 0, 'D': 1, 'L': 2, 'R': 3}[_edge_name_]

    @classmethod
    def _edge_index_to_name_dict_(cls):
        return {0: 'U', 1: 'D', 2: 'L', 3: 'R'}

    @classmethod
    def _edge_index_to_name_(cls, _edge_index_):
        return {0: 'U', 1: 'D', 2: 'L', 3: 'R'}[_edge_index_]

    @classmethod
    def _edge_pairing_(cls):
        return {'U': 'D', 'D': 'U', 'L': 'R', 'R': 'L'}

    @classmethod
    def _edge_index_pairing_(cls):
        return {0: 1, 1: 0, 2: 3, 3: 2}

    @classmethod
    def _corner_name_to_index_dict_(cls):
        return {'UL': 0, 'DL': 1, 'UR': 2, 'DR': 3,
                'LU': 0, 'LD': 1, 'RU': 2, 'RD': 3}

    @classmethod
    def _corner_name_to_index_(cls, _corner_name_):
        return {'UL': 0, 'DL': 1, 'UR': 2, 'DR': 3,
                'LU': 0, 'LD': 1, 'RU': 2, 'RD': 3}[_corner_name_]

    @classmethod
    def _corner_index_to_name_(cls, _corner_index_):
        return {0: 'UL', 1: 'DL', 2: 'UR', 3: 'DR'}[_corner_index_]

    @classmethod
    def _edge_corner_local_numbering_(cls, _edge_index_):
        """
        Notice the values are always increasing.

        """
        return {0: (0, 2),  # the Upper-side has corners locally numbered 0, 2.
                1: (1, 3),  # the Down-side has corners locally numbered 1, 3.
                2: (0, 1),  # the Left-side has corners locally numbered 0, 1.
                3: (2, 3)}[_edge_index_]

    @classmethod
    def _axis_indix_dict_(cls):
        return {'x': 0, 'y': 1}

    @classmethod
    def _edge_axis_distribution_(cls):
        """
        Here 'U':(0,0) means the U edge is perpendicular to the 0-axis and at the
        starting side. 'D':(0,-1) means the 'D' edge is perpendicular to the 0-axis but
        at the end side, so -1. And 'R':(1,-1) means the 'R' edge is perpendicular to
        1-axis and at the end side. And so on.

        """
        return {'U': (0, 0), 'D': (0, -1), 'L': (1, 0), 'R': (1, -1)}


