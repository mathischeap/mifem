


def parse_boundary_position_indicator(element):
    """

    Parameters
    ----------
    element :
        The node-element instance.

    Returns
    -------
    indicator :

        'corner|6NEF' : For example, the 8 corners of the crazy mesh.
        'corner-edge|2NEF.6NEB' : For example, at the 12 domain edge (except the domain corner) of
            the crazy mesh.
'
    """

    positions = element.positions

    boundary_positions = list()
    non_boundary_positions = list()

    for pos in positions:
        if pos[0].isnumeric():
            non_boundary_positions.append(pos)
        else:
            boundary_positions.append(pos)

    if len(non_boundary_positions) == 1: # must be at the boundary corner of domain corner mesh-element
        # for example, the 8 corner node-elements of the crazy mesh.
        return 'corner|' + non_boundary_positions[0]

    elif len(non_boundary_positions) == 2: # must be at the boundary corner-edge
        # for example, the node-elements at the 12 domain edges of the crazy mesh.
        return 'corner-edge|' + non_boundary_positions[0] + '.' + non_boundary_positions[1]

    elif len(non_boundary_positions) == 4: # may be one of multiple cases.
        # the simplest situation is at a surface.

        set0 = set(non_boundary_positions[0][-3:])
        set1 = set(non_boundary_positions[1][-3:])
        set2 = set(non_boundary_positions[2][-3:])
        set3 = set(non_boundary_positions[3][-3:])

        SHARE = set0 & set1 & set2 & set3

        if len(SHARE) == 1: # this node-element is at the middle of a boundary-surface.

            norm_direction = list(SHARE)[0]

            indicator = 'surface-middle|' + norm_direction
            for pos in non_boundary_positions:
                indicator += '-' + pos

            return indicator
        else:
            raise NotImplementedError()

    else:
        raise NotImplementedError()