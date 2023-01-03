from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.lib.distances as mdadist
import MDAnalysis.analysis.rms as rms

"""
DESCRIPTION
    This Python file contains the calculators that are called from measurements.py.
    Each function needs to input the selections and the options.
"""

def build_plane(
    positions,  # nx3, being n the number of atoms that build the plane
):
    """
    DESCRIPTION
        Function for obtaining the equation of the plane from a given list containing the positions of three atoms.
        It is used in the 'planar_ange' function.
    """
    from numpy import cross, dot

    if len(positions) == 3:
        v1 = positions[2] - positions[0]
        v2 = positions[1] - positions[0]

        a, b, c = cross(v1, v2)
        d = dot([a, b, c], positions[2])

        #        return array(a, b, c, d)
        return [a, b, c, d]


def planar_angle(
    plane_A,  # nx3 np array, being n the number of atoms for building the plane (min 3)
    plane_B,  # nx3 np array, being n the number of atoms for building the plane (min 3)
    units,  # [ deg | rad ]
    domain,  # [ 360 | 180 ] 360: (0, 360); 180: (-180, 180)
):
    from numpy import arccos, sqrt

    plane_A = build_plane(plane_A)
    plane_B = build_plane(plane_B)

    ang = arccos(
        (plane_A[0] * plane_B[0] + plane_A[1] * plane_B[1] + plane_A[2] * plane_B[2])
        / (
            sqrt(plane_A[0] ** 2 + plane_A[1] ** 2 + plane_A[2] ** 2)
            * sqrt(plane_B[0] ** 2 + plane_B[1] ** 2 + plane_B[2] ** 2)
        )
    )

    if units in ("rad", "radian", "radians"):

        if domain in (180, "180", "pi"):
            pass

        elif domain in (360, "360", "2pi"):
            from math import pi

            ang = (ang + pi) % pi

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        if domain in (180, "180", "pi"):
            ang = rad2deg(ang)

        elif domain in (360, "360", "2pi"):
            ang = (rad2deg(ang) + 360) % 360

    return ang


def distance(
    sel1, # selection containing one or more atoms
    sel2, # selection containing one or more atoms
    type # [min | max | com | cog ]
):
    """
    DESCRIPTION
        Function that takes two selections and gives:
            - minimum distance
            - maximum distance
            - distance between Centers of Mass (com)
            - distance between Centers of Geometry (cog)
    """
    from numpy import min as npmin
    from numpy import max as npmax

    if type == 'min':
        d = npmin(mdadist.distance_array(sel1, sel2, backend='OpenMP'))

    elif type == 'max':
        d = npmax(mdadist.distance_array(sel1, sel2, backend='OpenMP'))

    elif type == 'com':
        d = mdadist.distance_array(sel1.center_of_mass(), sel2.center_of_mass(), backend='OpenMP')

    elif type == 'cog':
        d = mdadist.distance_array(sel1.center_of_geometry(), sel2.center_of_geometry(), backend='OpenMP')

    return d


def dihedral(
    sel1,
    sel2,
    sel3,
    sel4,
    units, # [ rad | deg ]
    domain # [ 180 | 360 ]
):
    """
    DESCRIPTION
        Function that calculates the dihedral angle between four atoms.
        Dihedral angle is calculated by calculating the planar angle of the planes builds by the sel1, sel2, sel3 and sel2, sel3, sel4.

    OPTIONS (as arguments)
        - units: radians (rad) or degrees (deg)
        - domain: -180 to 180 º or 0 to 360º
    """

    d = mdadist.calc_dihedrals(sel1, sel2, sel3, sel4, backend="OpenMP")

    if units in ("rad", "radian", "radians", "pi"):
        from math import pi

        if domain in (360, "360", "2pi"):
            d = ( d + pi ) % pi
            return d
        elif domain in (180, "180", "pi"):
            return d

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        d = rad2deg(d)
        if domain in (360, "360", "2pi"):
            d = ( d + 360 ) % 360
            return d
        elif domain in (180, "180", "pi"):
            return d



def angle(
    sel1,
    sel2,
    sel3,
    units, # [ rad | deg ]
    domain # [ 180 | 360 ]
):
    """
    DESCRIPTION
        Function that calculates the angle between three atoms.
        Order of selection is important (sel2 is the vertex)

    OPTIONS (as arguments)
        - units: radians (rad) or degrees (deg)
        - domain: -180 to 180 º or 0 to 360º
    """

    a = mdadist.calc_angles(sel1, sel2, sel3, backend="OpenMP")

    if units in ("rad", "radian", "radians", "pi"):
        from math import pi

        if domain in (360, "360", "2pi"):
            a = ( a + pi ) % pi
            return a
        elif domain in (180, "180", "pi"):
            return a

    elif units in ("deg", "degree", "degrees"):
        from numpy import rad2deg

        a = rad2deg(a)
        if domain in (360, "360", "2pi"):
            a = ( a + 360 ) % 360
            return a
        elif domain in (180, "180", "pi"):
            return a