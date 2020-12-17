from ..exceptions import NotSingleAtomSelectionError

from MDAnalysis import Universe
import MDAnalysis.lib.distances as mdadist

# BOOLEAN CHECKERS
def dist_bool_output(dist, dist1, dist2=0, mode='lim'):
    """
    This function takes an input distance and the upper and lower limits or the central value and the tolerance and outputs
    if it satisfies or not the given criteria.

    modes:
        lim -> limit mode, it requires a max and a min value
        tol -> tolerance mode, it requires a central value and a min and max value
    """

    if mode == 'lim':
        if dist1 >= dist2:
            d_max, d_min = dist1, dist2
        elif dist2 > dist1:
            d_max, d_min = dist2, dist1
    elif mode == 'tol':
        d_max = dist1 + dist2
        d_min = dist1 - dist2

    return (dist <= d_max and dist > d_min)

def angle_bool_output(ang, ang1, ang2, mode='tol'):
    """
    This function takes an input angle and the upper and lower limits or the central value and the tolerance and outputs
    if it satisfies or not the given criteria.
    All the input values are transformed into the (0, 360) degrees format.
    The input values have to be in degrees, not in radians. They can be transformed using the np.rad2dreg() method, for example.

    modes:
        lim -> limit mode, it requires a max and a min value
        tol -> tolerance mode, it requires a central value and a min and max value
    """

    ang, ang1, ang2 = list(map(lambda ang: ((ang + 360) % 360), [ang, ang1, ang2]))

    if mode == 'lim':
        if ang1 >= ang2:
            a_max, a_min = ang1, ang2
        elif ang2 > ang1:
            a_max, a_min = ang2, ang1
    elif mode == 'tol':
        a_max = ang1 + ang2
        a_min = ang1 - ang2

    return (ang <= a_max and ang > a_min)



# DISTANCE, ANGLE AND DIHEDRAL ANGLES MEASURERS
def dist_measure(sel1, sel2):
    """
        This function outputs the minimum measured distance between the two input selections.
    """
    from numpy import min as npmin

    return npmin(mdadist.distance_array(sel1.positions, sel2.positions, backend='OpenMP'))

def dihe_measure(sel1, sel2, sel3, sel4):
    """
    This functions measures the dihedral angle between 4 specified atoms and returns the dihedral value between 0 and 360 degrees.
    The input selections have to be single atoms.
    """
    from numpy import rad2deg

    for sel in (sel1, sel2, sel3, sel4):
        if len(sel) != 1:
            raise NotSingleAtomSelectionError

    return ((float(rad2deg(mdadist.calc_dihedrals(sel1.positions, sel2.positions, sel3.positions, sel4.positions, backend='OpenMP')))) + 360) % 360

def ang_measure(sel1, sel2, sel3):
    """
    This functions measures the angle between 3 specified atoms and returns the value between 0 and 360 degrees.
    The input selections have to be single atoms.
    """
    from numpy import rad2deg

    for sel in (sel1, sel2, sel3):
        if len(sel) != 1:
            raise NotSingleAtomSelectionError

    return ((float(rad2deg(mdadist.calc_angles(sel1.positions, sel2.positions, sel3.positions, backend='OpenMP')))) + 360) % 360


######################

def dist_WATbridge(u, sel1, sel2, sel1_env=3, sel2_env=3):
    """
        This function takes two selection and looks for the nearest bridgin water between them.

        It requires the universe (it needs to select the environment of each selection), the two selections,
        the radius around each selection to look for waters (or the selections of the environment, which requires
        the updating=True argument).
        If no water has been found it will output 'None, None, None'.
    """
    from numpy import min as npmin

    if isinstance(sel1_env, (int, float)):
        res_list = list(sel1.residues)
        for i in range(len(res_list)):
            res_list[i] = str(res_list[i])[(str(res_list[i]).find(', ') + 2):(str(res_list[i]).find('>'))]

        sel1_env = list(u.select_atoms('around ' + str(sel1_env) +
                                            ' resid ' + " or resid ".join(res_list)).residues)

    elif isinstance(sel1_env, Universe):
        pass

    if isinstance(sel2_env, (int, float)):
        res_list = list(sel2.residues)
        for i in range(len(res_list)):
            res_list[i] = str(res_list[i])[(str(res_list[i]).find(', ') + 2):(str(res_list[i]).find('>'))]

        sel2_env = list(u.select_atoms('around ' + str(sel2_env) +
                                            ' resid ' + " or resid ".join(res_list)).residues)

    elif isinstance(sel2_env, Universe):
        pass

    WAT_list = []

    WAT_list = list(filter(lambda r: 'WAT' in str(r), list(set(sel1_env) & set(sel2_env))))
    WAT_list = str(WAT_list)[1:-1].replace('<Residue WAT, ', '')
    WAT_list = str(WAT_list).replace('>', '')
    WAT_list = WAT_list.replace(', ', ' ').split()

    if WAT_list != []:
        for WAT in WAT_list:
            dist1  = npmin(mdadist.distance_array(
                                (u.select_atoms('resid ' + WAT).positions),
                                sel1.positions,
                                backend='OpenMP'))
            dist2 = npmin(mdadist.distance_array(
                                u.select_atoms('resid ' + WAT).positions,
                                sel2.positions,
                                backend='OpenMP'))

            try :
                if dist1 < dist1_min and dist2 < dist2_min:
                    dist1_min  = dist1
                    dist2_min  = dist2
                    bridge_WAT = WAT

            except NameError:
                dist1_min  = dist1
                dist2_min  = dist2
                bridge_WAT = WAT

        return dist1_min, dist2_min, bridge_WAT

    else :
        return None, None, None

def dist_plane(sel, plane1, plane2=None, plane3=None):
    """
        This functions takes one input atom and 3 other input atoms that constitute a plane
        and outputs the distance value and sign.
        #    There are two sets of inputs: the sel, which is the atom that will be measured and which can be input as
        #    a single atom selection, as a selection of multiple atoms or as a list of multiple selections; and the
        #    plane1 (plane2, plane3), which are the three atoms that
        Input selection has to be a selection of the atoms to be measured, while plane atoms can be either a
        selection of 3 atoms or 3 selection of 1 atom.
        The output is a float value or a list of float values depending on the length of the input selection
    """

    def plane_eq(plane1, plane2, plane3):
        from numpy import cross, dot

        v1 = plane2.position - plane1.position
        v2 = plane3.position - plane1.position

        a, b, c = cross(v1, v2)
        d = dot([a, b, c], plane3.position)

        return [a, b, c, d]

    def distance(point, eq):
        from math import sqrt

        return ((eq[0]*point[0] + eq[1]*point[1] + eq[2]*point[2] - eq[3])/(sqrt((eq[0] ** 2) + (eq[1] ** 2) + (eq[2] ** 2))))

    if len(sel) == 1:
        sel = sel.atoms[0]

    else :
        raise NotSingleAtomSelectionError


    if len(plane1) == 1 and len(plane2) == 1 and len(plane3) == 1:
        plane3 = plane3.atoms[0]
        plane2 = plane2.atoms[0]
        plane1 = plane1.atoms[0]

    elif len(plane1) == 2 and len(plane2) == 1 and plane3 == None:
        plane3 = plane2.atoms[0]
        plane2 = plane1.atoms[1]
        plane1 = plane1.atoms[0]

    elif len(plane1) == 1 and len(plane2) == 2 and plane3 == None:
        plane3 = plane2.atoms[1]
        plane2 = plane2.atoms[0]
        plane1 = plane1.atoms[0]

    elif len(plane1) == 3 and plane2 == None and plane3 == None:
        plane3 = plane1.atoms[2]
        plane2 = plane1.atoms[1]
        plane1 = plane1.atoms[0]

    else :
        raise NotSingleAtomSelectionError


    return distance(sel.position, plane_eq(plane1, plane2, plane3))
