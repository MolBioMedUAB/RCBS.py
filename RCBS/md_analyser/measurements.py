# from asyncio import selector_events
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.lib.distances as mdadist
from ..exceptions import (
    NotExistingInteraction,
    NotSingleAtomSelectionError,
    NotExistingInteraction,
    OutputFormatNotAvailable,
)
from .selections import selection
from numpy import min as npmin
from numpy import array


class Measurements:
    """
    DESCRIPTION:
        Class containing all the possible analysis of a trajectory that RCBS is able to carry out. It creates a dictionary of results (self.results) containing all the results of the configurated measurements (by add_* functions) and calculated using run_measure function, and creates a dictionary of boolean values (self.boolean) containing the result of applying a certain criterium to a certain measurement.

    FUNCTIONS:
        - add_distance: measure the distance between two atoms or two groups of atoms (and return the shorter one)
        - add_dihedral: measure the dihedral angle between four atoms
        - add_angle:    measure the angle between three atoms
        - add_contacts: measure the number of contacts at a given distance from a given selection of one or more atoms
        - run_measure:  run the all the measurements previously added
        - run_boolean:  check if the measured parameters satisfy the given criteria
        - config_saver: save the dictionary in JSON or YAML format containing all the measurements to carry out with the run_measure function

    INPUT:
        - u: preloaded MDAnalysis' universe

    OUTPUT:
        - self.results: dictionary containing results keyed by the given name
        - self.boolean: dictionary containing the boolean results keyed by the given name of the criteria applied to self.results
    """

    def __init__(self, u):
        self.measurements = []
        self.universe = u
        self.results = {}
        self.boolean = {}

    def add_distance(self, name, sel1, sel2, type="min"):
        """
        DESCRIPTION:
            This function outputs the minimum measured distance between the two input selections or coordinates or their combination.

        INPUT:
            - Name of the measurement
            - Two selections, which can contain more than one atom
            - type: min [default], com (center of mass) or cog (center of geometry)

        OUTPUT:
            - Shorter distance between sel1 and sel2 (in ang)
        """
        if type.lower() not in ("min", "com", "cog"):
            raise NotAvailableOptionError

        self.measurements.append(
            {
                "name": name,
                "type": "distance",
                "sel": [sel1, sel2],
                "options": {"type": type},
            }
        )

    def add_dihedral(self, name, sel1, sel2, sel3, sel4, units="degree", domain=360):
        """
        DESCRIPTION:
            This functions measures the dihedral angle between 4 specified atoms and returns the dihedral value between 0 and 360 degrees.
            The input selections have to be single atoms.

        OPTIONS:
            - units: option for selecting the output units of the dihedral
                - degree
                - rad

            - domain: option for specifying the domain of the output measures
                - 180, pi: option for -180,180 domain
                - 360, 2pi: option for 0,360 domain. Default option

        INPUT:
            - Name of the measurement
            - Selection of four atoms in four different AtomGroups. They have to be input with the correct order

        OUTPUT:
            - Dihedral angle between the input atoms
        """

        for sel in (sel1, sel2, sel3, sel4):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError

        units = units.lower()
        domain = str(domain).lower()

        if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
            units = "degree"

        if domain not in (180, "180", 360, "360", "pi", "2pi"):
            domain = "360"

        self.measurements.append(
            {
                "name": name,
                "type": "dihedral",
                "sel": [sel1, sel2, sel3, sel4],
                "options": {"units": units, "domain": domain},
            }
        )

    def add_angle(self, name, sel1, sel2, sel3, units="deg", domain=360):
        """
        DESCRIPTION:
            This functions measures the angle between 3 specified atoms and returns the value between 0 and 360 degrees.
            The input selections have to be single atoms.

        OPTIONS:
            - Name of the measurement
            - units: option for selecting the output units of the dihedral
                - degree
                - rad

            - domain: option for specifying the domain of the output measures
                - 180, pi: option for -180,180 domain
                - 360, 2pi: option for 0,360 domain. Default option

        INPUT:
            - Selection of four atoms in three different AtomGroups. They have to be input with the correct order

        OUTPUT:
            - Angle between the input atoms
        """

        for sel in (sel1, sel2, sel3):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError

        units = units.lower()
        domain = str(domain).lower()

        if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
            units = "degree"

        if domain not in (180, "180", 360, "360", "pi", "2pi"):
            domain = "360"

        self.measurements.append(
            {
                "name": name,
                "type": "angle",
                "sel": [sel1, sel2, sel3],
                "options": {"units": units, "domain": domain},
            }
        )

    def add_contacts(self, name, sel, sel_env=3, interactions="all", include_WAT=False):
        """
        DESCRIPTION:
            This function takes a Universe, a selection and a radius and returns the list of residues nearer than the specified radius.

        INPUT:
            - Name of the measurement
            - u            -> MDAnalysis Universe
            - sel          -> selection of central atoms. It has to be an AtomGroup
            - sel_env      -> radius (in ang)
            - interactions -> type of interactions to be considered (all, polar, nonpolar, donorHbond, none). Custom
                              interactions can be also analysed by passing a list of residues names

        OUTPUT:
            - List of dictionaries containing the name and number of all interacting residues
        """

        sel_env = self.universe.select_atoms(
            "around %s group select" % sel_env, select=sel, updating=True
        )

        if isinstance(interactions, str):
            if interactions not in ("all", "polar", "nonpolar", "donorHbond", "none"):
                raise NotExistingInteraction

            else:
                if interactions == "all":
                    interactions = [
                        "ARG",
                        "HIS",
                        "HID",
                        "HIE",
                        "HIP",
                        "LYS",
                        "ASP",
                        "ASH",
                        "GLU",
                        "GLH",
                        "SER",
                        "THR",
                        "ASN",
                        "GLN",
                        "CYS",
                        "SEC",
                        "GLY",
                        "PRO",
                        "ALA",
                        "ILE",
                        "LEU",
                        "MET",
                        "PHE",
                        "TRP",
                        "TYR",
                        "VAL",
                    ]
                elif interactions == "polar":
                    interactions = [
                        "ARG",
                        "HIS",
                        "HID",
                        "HIE",
                        "HIP",
                        "LYS",
                        "ASP",
                        "ASH",
                        "GLU",
                        "GLH",
                        "SER",
                        "THR",
                        "ASN",
                        "GLN",
                        "CYS",
                        "SEC",
                        "TYR",
                    ]

                elif interactions == "nonpolar":
                    interactions = [
                        "CYS",
                        "SEC",
                        "GLY",
                        "PRO",
                        "ALA",
                        "ILE",
                        "LEU",
                        "MET",
                        "PHE",
                        "TRP",
                        "TYR",
                        "VAL",
                    ]

                elif interactions == "donorHbond":
                    interactions = [
                        "ARG",
                        "HID",
                        "HIE",
                        "HIP",
                        "LYS",
                        "ASH",
                        "GLH",
                        "SER",
                        "THR",
                        "ASN",
                        "GLN",
                        "CYS",
                        "SEC",
                        "GLY",
                        "PRO",
                        "TYR",
                    ]

                elif interactions == "none":
                    interactions = []

        elif isinstance(interactions, list):

            pass

        if include_WAT == True:
            interactions += ["WAT"]

        self.measurements.append(
            {
                "name": name,
                "type": "contacts",
                "sel": [sel, sel_env],
                "options": {
                    "interactions": interactions,
                },
            }
        )

    def add_distWATbridge(self, name, sel1, sel2, sel1_env=3, sel2_env=3):
        """
        DESCRIPTION
           This function takes a Universe, two selections and the size of their environments and returns the nearest bridging water between the two selections and the distance to both of them.

        INPUT:
            - Name of the measurement
            - u            -> MDAnalysis Universe
            - sel1         -> selection of first set of central atoms. It has to be an AtomGroup
            - sel2         -> selection of second set of central atoms. It has to be an AtomGroup
            - sel1_env      -> radius around the first set of central atoms (in ang)
            - sel2_env      -> radius around the first set of central atoms (in ang)

        OUTPUT:
            - List of dictionaries containing the number of the bridging water and the smallest distance to each of the selection sets.
        """

        sel1_env = self.universe.select_atoms(
            "resname WAT and around %s group select" % sel1_env,
            select=sel1,
            updating=True,
        )

        sel2_env = self.universe.select_atoms(
            "resname WAT and around %s group select" % sel2_env,
            select=sel2,
            updating=True,
        )

        self.measurements.append(
            {
                "name": name,
                "type": "distWATbridge",
                "sel": [sel1, sel2, sel1_env, sel2_env],
                "options": None,
            }
        )

    def config_saver(self, config_filename, verbose=True):
        """
        DESCRIPTION:
            Function for saving the configuration of a measurement run into a file (in json or yaml format).

        INPUT:
            - self (containing the measurements list of dictionaries) and the file name

        OUTPUT:
            - the dictionary saved in the file
        """

        if config_filename.split(".")[-1].lower() not in ("json", "jsn", "yaml", "yml"):
            raise OutputFormatNotAvailable

        config = []
        for l in self.measurements:
            config.append(l.copy())

        for measurement in config:

            sels = []
            for s in measurement["sel"]:
                if isinstance(s, AtomGroup):
                    sels.append(s.indices)

                    for m in range(len(sels)):
                        sels[m] = list(sels[m])
                        for n in range(len(sels[m])):
                            sels[m][n] = int(sels[m][n] + 1)

                measurement["sel"] = sels

            measurement["sel_type"] = "at_num"

        if config_filename.split(".")[-1].lower() in ("json", "jsn"):
            from json import dump

            f = open(config_filename, "w")
            dump(config, f)
            f.close()

            if verbose == True:
                print("Configuration saved in", config_filename)

            del dump

        elif config_filename.split(".")[-1].lower() in ("yaml", "yml"):
            from yaml import dump

            f = open(config_filename, "w")
            dump(config, f)
            f.close()

            if verbose == True:
                print("Configuration saved in", config_filename)

            del dump

        return config

    def run_measure(self, save_output=False, verbose=True):
        """
        DESCRIPTION:
            Function for runninng all the configured measurments on a given trajectory (loaded as self.universe). It can take also a configuration stored in a file instead of taking the in-situ configurated measurement.

        OPTIONS:
            - save_output: Pseudoboolean value for storing the resutls as a file. False means no storing, any other string with the json or yaml extension means save the results in a file called as given.

        INPUT:
            - self: containing self.universe and self.measurements (if not loading a configuration file)

        OUPTUT:
            - self.results: dictionary containing arrays of the results keyed by the given name of the configuration.
        """

        for measurement in self.measurements:
            self.results[measurement["name"]] = []

        for ts in self.universe.trajectory:

            for measurement in self.measurements:

                if measurement["type"] == "distance":
                    if measurement["options"]["type"] == "min":
                        self.results[measurement["name"]].append(
                            npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][0].positions),
                                    array(measurement["sel"][1].positions),
                                    backend="OpenMP",
                                )
                            )
                        )

                    elif measurement["options"]["type"] == "com":
                        self.results[measurement["name"]].append(
                            npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][0].center_of_mass()),
                                    array(measurement["sel"][1].center_of_mass()),
                                    backend="OpenMP",
                                )
                            )
                        )

                    elif measurement["options"]["type"] == "cog":
                        self.results[measurement["name"]].append(
                            npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][0].center_of_geometry()),
                                    array(measurement["sel"][1].center_of_geometry()),
                                    backend="OpenMP",
                                )
                            )
                        )

                elif measurement["type"] == "dihedral":

                    if measurement["options"]["units"] in ("rad", "radian", "radians"):
                        if measurement["options"]["domain"] in (180, "180", "pi"):
                            self.results[measurement["name"]].append(
                                float(
                                    mdadist.calc_dihedrals(
                                        measurement["sel"][0].positions,
                                        measurement["sel"][1].positions,
                                        measurement["sel"][2].positions,
                                        measurement["sel"][3].positions,
                                        backend="OpenMP",
                                    )
                                )
                            )
                        elif measurement["options"]["domain"] in (360, "360", "2pi"):
                            from math import pi

                            self.results[measurement["name"]].append(
                                (
                                    float(
                                        mdadist.calc_dihedrals(
                                            measurement["sel"][0].positions,
                                            measurement["sel"][1].positions,
                                            measurement["sel"][2].positions,
                                            measurement["sel"][3].positions,
                                            backend="OpenMP",
                                        )
                                    )
                                    + pi
                                )
                                % pi
                            )

                    elif measurement["options"]["units"] in (
                        "deg",
                        "degree",
                        "degrees",
                    ):
                        from numpy import rad2deg

                        if measurement["options"]["domain"] in (180, "180", "pi"):
                            self.results[measurement["name"]].append(
                                float(
                                    rad2deg(
                                        mdadist.calc_dihedrals(
                                            measurement["sel"][0].positions,
                                            measurement["sel"][1].positions,
                                            measurement["sel"][2].positions,
                                            measurement["sel"][3].positions,
                                            backend="OpenMP",
                                        )
                                    )
                                )
                            )
                        elif measurement["options"]["domain"] in (360, "360", "2pi"):
                            from math import pi

                            self.results[measurement["name"]].append(
                                (
                                    float(
                                        rad2deg(
                                            mdadist.calc_dihedrals(
                                                measurement["sel"][0].positions,
                                                measurement["sel"][1].positions,
                                                measurement["sel"][2].positions,
                                                measurement["sel"][3].positions,
                                                backend="OpenMP",
                                            )
                                        )
                                    )
                                    + 360
                                )
                                % 360
                            )

                elif measurement["type"] == "angle":

                    if measurement["options"]["units"] in ("rad", "radian", "radians"):
                        if measurement["options"]["domain"] in (180, "180", "pi"):
                            self.results[measurement["name"]].append(
                                float(
                                    mdadist.calc_angles(
                                        measurement["sel"][0].positions,
                                        measurement["sel"][1].positions,
                                        measurement["sel"][2].positions,
                                        backend="OpenMP",
                                    )
                                )
                            )
                        elif measurement["options"]["domain"] in (360, "360", "2pi"):
                            from math import pi

                            self.results[measurement["name"]].append(
                                (
                                    float(
                                        mdadist.calc_angles(
                                            measurement["sel"][0].positions,
                                            measurement["sel"][1].positions,
                                            measurement["sel"][2].positions,
                                            backend="OpenMP",
                                        )
                                    )
                                    + pi
                                )
                                % pi
                            )

                    elif measurement["options"]["units"] in (
                        "deg",
                        "degree",
                        "degrees",
                    ):
                        from numpy import rad2deg

                        if measurement["options"]["domain"] in (180, "180", "pi"):
                            self.results[measurement["name"]].append(
                                float(
                                    rad2deg(
                                        mdadist.calc_angles(
                                            measurement["sel"][0].positions,
                                            measurement["sel"][1].positions,
                                            measurement["sel"][2].positions,
                                            backend="OpenMP",
                                        )
                                    )
                                )
                            )
                        elif measurement["options"]["domain"] in (360, "360", "2pi"):
                            from math import pi

                            self.results[measurement["name"]].append(
                                (
                                    float(
                                        rad2deg(
                                            mdadist.calc_angles(
                                                measurement["sel"][0].positions,
                                                measurement["sel"][1].positions,
                                                measurement["sel"][2].positions,
                                                backend="OpenMP",
                                            )
                                        )
                                    )
                                    + 360
                                )
                                % 360
                            )

                elif measurement["type"] == "contacts":
                    names = list(measurement["sel"][1].resnames)
                    ids = list(measurement["sel"][1].resids)

                    dict_ = {}

                    for i in range(len(ids)):
                        if ids[i] not in measurement["sel"][0].residues.resindices:
                            if str(names[i]) in measurement["options"]["interactions"]:
                                dict_[int(ids[i])] = names[i]

                    #                    for i, n in zip(ids, names):
                    #                        if i not in measurement["sel"][1].residues.resindices:
                    #                            if str(n)[:3] in measurement["options"]["interactions"]:
                    #                                dict_[int(str(i))] = n

                    self.results[measurement["name"]].append(dict_)

                elif measurement["type"] == "distWATbridge":
                    if (
                        "WAT" not in measurement["sel"][2].resnames
                        or "WAT" not in measurement["sel"][3].resnames
                    ):  # No WAT in one or both selections' environment
                        self.results[measurement["name"]].append([None, None, None])

                    elif (
                        len(
                            set(measurement["sel"][2].resids)
                            & set(measurement["sel"][3].resids)
                        )
                        == 0
                    ):  # No WAT present in both environments
                        self.results[measurement["name"]].append([None, None, None])

                    else:  # WAT residues in both environments and at least one of them coincident
                        WATs = list(
                            set(measurement["sel"][2].resids)
                            & set(measurement["sel"][3].resids)
                        )

                        for wat in WATs:
                            selWAT = selection(
                                self.universe, int(wat), sel_type="res_num"
                            )

                            dist1_ = npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][0].positions),
                                    array(selWAT.positions),
                                    backend="OpenMP",
                                )
                            )

                            dist2_ = npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][1].positions),
                                    array(selWAT.positions),
                                    backend="OpenMP",
                                )
                            )

                            try:
                                if (dist1_ + dist2_) / 2 < (
                                    dist1 + dist2
                                ) / 2:  # Closest WAT is the one with the shortest average distance to both of the sels
                                    dist1, dist2 = dist1_, dist2_
                                    closestWAT = wat

                                else:
                                    pass

                            except NameError:
                                dist1, dist2 = dist1_, dist2_
                                closestWAT = wat

                    self.results[measurement["name"]].append([closestWAT, dist1, dist2])

                    del selWAT

        if save_output != False:

            if save_output.split(".")[-1] in ("json", "jsn"):
                from json import dump

                f = open(save_output, "w")
                dump(self.results, f)
                f.close()
                if verbose == True:
                    print("Results saved in", save_output)

            elif save_output.split(".")[-1] in ("yaml", "yml"):
                from yaml import dump

                f = open(save_output, "w")
                dump(self.results, f)
                f.close()

                if verbose == True:
                    print("Results saved in", save_output)

        return

    def run_boolean(self, *bool_configs, combine=True, save_output=False, verbose=True):
        """
        DESCRIPTION:
            Function for checking if the measured results satisfy the given criteria.

        CHECKERS:
            List of con
            {
                measure_name : , #mandatory. Has to be the same from the results dictionary.
                measure_type : dist | ang | dihe, # default is dist
                mode : lim | tol, #default is lim
                ref_val1 : , #mandatory
                ref_val2 :  # default is 0
            }
            EX: {measure_name : , measure_type : dist, mode : lim, ref_val1

        INPUT:
            - self: containing the self.results dictionary
            - bool_configs: list of dictionaries containing the information about the criteria to apply. The 'CHECKERS' section explains the structure of the dictionaries.

        OPTIONS:
            - combine: boolean argument. If true it adds an array containing the boolean per-frame values of the combination of all the other criteria.
                    The combination is calculated using the and logivcal operation.
            - save_output: Pseudoboolean value for storing the resutls as a file. False means no storing, any other string with the json or yaml extension
                    means save the results in a file called as given.
        """

        def bool_configs_checker(bool_configs):
            """
            DESCRIPTION
                Function for checking if the dictionaries containing the information about the criteria fit contain all the required information.

            INPUT:
                - List of dictionaries
            """

            counter = 0
            for c in bool_configs:

                if "measure_name" not in list(c.keys()):
                    c["measure_name"] = "boolean_" + str(counter)

                if "measure_type" not in list(c.keys()):
                    c["measure_type"] = "dist"
                else:
                    if c["measure_type"].lower() not in ("dist", "ang", "dihe"):
                        c["measure_type"] = "dist"

                if "mode" not in list(c.keys()):
                    c["mode"] = "lim"
                else:
                    if c["mode"].lower() not in ("lim", "tol"):
                        c["mode"] = "lim"

                if "ref_val1" not in list(c.keys()):
                    while True:
                        try:
                            c["ref_val1"] = float(
                                input(
                                    "Input max (if lim mode) or central (if tol mode) value [float] for %s: "
                                    % c["measure_name"]
                                )
                            )
                            break

                        except ValueError:
                            print("Input a number (int or float).")
                            continue

                if "ref_val2" not in list(c.keys()):
                    if c["mode"] == "lim":
                        c["ref_val2"] = 0

                    elif c["mode"] == "tol":
                        while True:
                            try:
                                c["ref_val2"] = float(
                                    input(
                                        "Input tolerance value (val_ref2) (if tol mode) value [float] for %s: "
                                        % c["measure_name"]
                                    )
                                )
                                break

                            except ValueError:
                                print("Input a number (int or float).")
                                continue

                counter += 1

            return bool_configs

        bool_configs = bool_configs_checker(bool_configs)

        for b in bool_configs:

            if b["mode"] == "lim":
                if b["ref_val1"] >= b["ref_val2"]:
                    max_val = b["ref_val1"]
                    min_val = b["ref_val2"]
                elif b["ref_val1"] < b["ref_val2"]:
                    max_val = b["ref_val2"]
                    min_val = b["ref_val1"]

                self.boolean[b["measure_name"]] = []
                for r in self.results[b["measure_name"]]:
                    self.boolean[b["measure_name"]].append(
                        bool(r <= max_val and r > min_val)
                    )

            elif b["mode"] == "tol":
                max_val = b["ref_val1"] + b["ref_val2"]
                min_val = b["ref_val1"] - b["ref_val2"]

                if b["measure_type"] == "dist":
                    self.boolean[b["measure_name"]] = []
                    for r in self.results[b["measure_name"]]:
                        self.boolean[b["measure_name"]].append(
                            bool(r <= max_val and r > min_val)
                        )

                elif b["measure_type"] in ("ang", "dihe"):
                    # max_val, min_val = list(map(lambda ang: ((ang + 360) % 360), [max_val, min_val]))

                    if max_val > 360 and min_val < 360:
                        self.boolean[b["measure_name"]] = []
                        for r in self.results[b["measure_name"]]:
                            r_ = (r + 360) % 360
                            self.boolean[b["measure_name"]].append(
                                bool(
                                    (r_ > min_val and r_ < 360)
                                    or (r_ > 0 and r_ < (max_val - 360))
                                )
                            )

                    elif min_val < 0 and max_val > 0:
                        self.boolean[b["measure_name"]] = []
                        for r in self.results[b["measure_name"]]:
                            r_ = (r + 360) % 360
                            self.boolean[b["measure_name"]].append(
                                bool(
                                    (r_ > 0 and r_ < max_val)
                                    or (r_ > (360 + min_val) and r_ < 360)
                                )
                            )

                    elif (min_val < 0 and max_val < 0) or (
                        min_val > 360 and max_val > 360
                    ):

                        max_val, min_val = list(
                            map(lambda ang: ((ang + 360) % 360), [max_val, min_val])
                        )

                        self.boolean[b["measure_name"]] = []
                        for r in self.results[b["measure_name"]]:
                            r_ = (r + 360) % 360
                            self.boolean[b["measure_name"]].append(
                                bool(r_ <= max_val and r_ > min_val)
                            )

        if combine == True:
            self.boolean["combination"] = []
            for f in range(len(self.boolean[list(self.boolean.keys())[0]])):

                boolean = True
                for key in list(self.boolean.keys())[:-1]:
                    boolean = bool(boolean and self.boolean[key][f])

                self.boolean["combination"].append(boolean)

        if save_output != False:

            if save_output.split(".")[-1] in ("json", "jsn"):
                from json import dump

                f = open(save_output, "w")
                dump(self.boolean, f)
                f.close()
                if verbose == True:
                    print("Boolean results saved in", save_output)

            elif save_output.split(".")[-1] in ("yaml", "yml"):
                from yaml import dump

                f = open(save_output, "w")
                dump(self.boolean, f)
                f.close()

                if verbose == True:
                    print("Boolean results saved in", save_output)

        return
