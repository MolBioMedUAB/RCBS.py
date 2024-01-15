# from asyncio import selector_events
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.lib.distances as mdadist
import MDAnalysis.analysis.rms as rms

from tqdm.autonotebook import tqdm

from numpy import min as npmin
from numpy import max as npmax
from numpy import array

from ..exceptions import (
    NotExistingInteraction,
    NotSingleAtomSelectionError,
    NotExistingInteraction,
    OutputFormatNotAvailable,
    NotAvailableOptionError,
    NotThreeAtomsSelectionError,
    EmptyMeasurementsError,
)
from .selections import selection
from .calculators import *


class Measurements:
    """
    DESCRIPTION:
        Class containing all the possible analysis of a trajectory that RCBS is able to carry out. It creates a dictionary of results (self.results) containing all the results of the configurated measurements (by add_* functions) and calculated using run_measure function, and creates a dictionary of boolean values (self.boolean) containing the result of applying a certain criterium to a certain measurement.

    FUNCTIONS:
        - add_distance:     measure the distance between two atoms or two groups of atoms (and return the shorter one)
        - add_dihedral:     measure the dihedral angle between four atoms
        - add_angle:        measure the angle between three atoms
        - add_contacts:     measure the number of contacts at a given distance from a given selection of one or more atoms or the residue-wise contacts of all the protein
        - add_RMSD:         measure the RMSD for a selection of each frame against a reference (a given structure or the first frame of the trajectory)
        - add_planar_angle: measure the angle between two planes defined by three atoms each one.
        - add_pKa:          measure the residue-wise pKa.
        - run_measure:      run the all the measurements previously added
        - run_boolean:      check if the measured parameters satisfy the given criteria
        - config_saver:     save the dictionary in JSON or YAML format containing all the measurements to carry out with the run_measure function

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
            - Shorter distance between sel1 and sel2 (in ang) or distances between COMs or COGs.
        """
        if type.lower() not in ("min", "max", "com", "cog"):
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

    def add_planar_angle(self, name, sel1, sel2, units="deg", domain=360):
        """
        DESCRIPTION:
            This function measures the angle between two planes specified by three atoms each one and returns the angle.
            The input selections have to contain three atoms.

        OPTIONS:
            - Name of the measurement
            - units: option for selecting the output units of the dihedral
                - degree
                - rad

            - domain: option for specifying the domain of the output measures
                - 180, pi: option for -180,180 domain
                - 360, 2pi: option for 0,360 domain. Default option

        INPUT:
            - Selection of two sets of three atoms in two different AtomGroups.

        OUTPUT:
            - Angle between the input atoms
        """

        for sel in (sel1, sel2):
            if len(sel) != 3:
                raise NotThreeAtomsSelectionError

        units = units.lower()
        domain = str(domain).lower()

        if units not in ("deg", "degree", "degrees", "rad", "radian", "radians"):
            units = "degree"

        if domain not in (180, "180", 360, "360", "pi", "2pi"):
            domain = "360"

        self.measurements.append(
            {
                "name": name,
                "type": "planar_angle",
                "sel": [sel1, sel2],
                "options": {"units": units, "domain": domain},
            }
        )

    def add_contacts(self, name, sel, sel_env=3, interactions="all", include_WAT=False, out_format='new', measure_distances=True):
        """
        DESCRIPTION:
            This function takes a Universe, a selection and a radius and returns the list of residues nearer than the specified radius.

        INPUT:
            - Name of the measurement
            - sel          -> selection of central atoms. It has to be an AtomGroup (only option for old engine) or the 'protein' string for 
                              measuring distances between all contacting residues
            - sel_env      -> radius (in ang)
            - interactions -> type of interactions to be considered (all, polar, nonpolar, donorHbond, none). Custom
                              interactions can be also analysed by passing a list of residues names
            - out_format   -> [ '0.3'/'new'/'n' | '0.2'/'old'/'o' ] Format of the output. 'old' corresponds to the old logics (versions 0.0 to 0.2)
                              and is kept for compatibility reasons.

        OUTPUT:
            - List of dictionaries containing the name and number of all interacting residues
        """

        #sel_env = self.universe.select_atoms(
        #    "around %s group select" % sel_env, select=sel, updating=True
        #)

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

        else :
            raise NotExistingInteraction
        

        if include_WAT == True:
            interactions += ["WAT", "HOH"]

        if isinstance(sel, AtomGroup):
            mode = 'selection'

            sel_env = self.universe.select_atoms(
                "around %s group select" % str(sel_env), select=sel, updating=True
            )

            if str(out_format).lower() in ['0.2', 'old', 'o']:
                out_format = 'old'
                if measure_distances:
                    print("Distances can not been calculated using the old output format. \
                          'measure_distances' has been set to False.")
                    out_format = 'new'
                    
            elif str(out_format).lower() in ['0.3', 'new', 'n']:
                out_format = 'new'

            else :
                print('The selected out format does not exist. The new format has been selected instead.')
                out_format = 'new'

            self.measurements.append(
            {
                "name": name,
                "type" : "contacts",
                "mode": mode,
                "sel": [sel, sel_env],
                "options": {
                    "interactions"  : interactions,
                    "measure_dists" : measure_distances,
                    "out_format"    : out_format
                },
            }
            )

        if isinstance(sel, str) and sel.lower() == 'protein':
            print('hey')
            
            mode = 'protein'
            #sel = self.universe.select_atoms('protein')
            sel = self.universe.select_atoms('all')

            self.measurements.append(
            {
                "name": name,
                "type" : "contacts",
                "mode": mode,
                "sel": [sel, sel_env],
                "options": {
                    "interactions"  : interactions,
                    "measure_dists" : measure_distances,
                    "out_format"    : 'new'
                },
            }
            )

    def add_RMSD(self, name, sel, ref=None, superposition=True):
        """
        DESCRIPTION:
            This function outputs the RMSD of a selection

        INPUT:
            - Name of the measurement
            - Selection
            - ref: selection of the reference universe. If not provided, the first frame will be used as the reference.
            - superposition [bool]: compute the RMSD of aligned

        OUTPUT:
            - Array of RMSDs of each frame against a reference
        """

        if isinstance(ref, type(None)):
            self.universe.trajectory[0]
            ref = sel.positions - sel.center_of_mass()

        elif isinstance(ref, AtomGroup):
            ref = ref.positions - ref.center_of_mass()

        self.measurements.append(
            {
                "name": name,
                "type": "rmsd",
                "sel": sel,
                "ref": ref,
                "options": {"superposition": superposition},
            }
        )

    def add_distWATbridge(self, name, sel1, sel2, sel1_rad=3, sel2_rad=3):
        """
        DESCRIPTION
           This function takes a Universe, two selections and the size of their environments and returns the nearest bridging water between the two selections and the distance to both of them.

        INPUT:
            - Name of the measurement
            - u            -> MDAnalysis Universe
            - sel1         -> selection of first set of central atoms. It has to be an AtomGroup
            - sel2         -> selection of second set of central atoms. It has to be an AtomGroup
            - sel1_rad      -> radius around the first set of central atoms (in ang)
            - sel2_rad      -> radius around the first set of central atoms (in ang)

        OUTPUT:
            - List of dictionaries containing the number of the bridging water and the smallest distance to each of the selection sets.
        """

        #sel1_rad, sel2_rad = sel1_env, sel2_env

        sel1_env = self.universe.select_atoms(
            "resname WAT and around %s group select" % sel1_rad,
            select=sel1,
            updating=True,
        )

        sel2_env = self.universe.select_atoms(
            "resname WAT and around %s group select" % sel2_rad,
            select=sel2,
            updating=True,
        )

        self.measurements.append(
            {
                "name": name,
                "type": "distWATbridge",
                "sel": [sel1, sel2, sel1_env, sel2_env, sel1_rad, sel2_rad],
                "options": None,
            }
        )

    def add_pKa(self, name, excluded_ions=["Na+", "Cl-"], pka_ref='neutral', pdb_folder='.pka', keep_pdb=False, keep_pka=False):
        """
        DESCRIPTION:
            This function allows the prediction of the pKa using PROpKa3 of the protein for each frame.

        INPUT:
            - name:             name of the measurement
            - excluded_ions:    list of solvent ions names that belong to solvent. Default are Na+ and Cl-.
            - pka_ref:          reference to calculate pKa. Default is neutral. [ neutral | low-pH]
            - keep_pdb:         trigger for keeping generated pdbs. Default is False. [ True | False ]
            - keep_pka:         trigger for keeping generated .pka file. Default is False. [ True | False ]

        OUPUT:
            - Per-frame array of dicts with shape { residue : pKa }
        """

        print('WARNING!!:')
        print('  Take into account that pKa calculation with PROpKA requires')
        print('  the generation of a PDB for each analysed frame.')
        print('  This can be very time consuming for long trajectory, so consider')
        print('  using the step option in run_measure() to reduce')
        print('  the amount of analysed structures.')

        excluded_ions = ' or resname '.join(excluded_ions)

        self.measurements.append(
            {
                "name" : name,
                "type" : "pka",
                "sel"  : self.universe.select_atoms("not (resname WAT or resname HOH or resname " + excluded_ions + ')'),
                "options" :
                    { 
                        "pka_ref" : pka_ref,
                        "pdb_folder" : pdb_folder,
                        "keep_pdb" : keep_pdb,
                        "keep_pka" : keep_pka
                    },
            }
        )
 
    def remove_measurement(self, name):
        if isinstance(name, str): name = [name] 

        for measurement_index in range(len(self.measurements)):
            for name_ in name:
                if self.measurements[measurement_index]['name'] == name_:
                    self.measurements.pop(measurement_index)

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

    def run_measure(self, step=1, first=1, last=-1, save_output=False, verbose=True, previous_measurements='overwrite'):
        """
        DESCRIPTION:
            Function for runninng all the configured measurments on a given trajectory (loaded as self.universe). It can take also a configuration stored in a file instead of taking the in-situ configurated measurement.

        OPTIONS:
            - save_output:              Pseudoboolean value for storing the resutls as a file. False means no storing, any other string with the json or yaml extension means save the results in a file called as given.
            - step:                     Frames to jump during the analysis. Default is 1, so all the trajectory will be analysed.
            - first:                    First frame to start the analysis. Default is 0.
            - last:                     Last frame to analyse (included). Default is last frame of trajectory
            - previous_measurements:    [ 'overwrite'/'o' | 'append'/'a' ] Way of dealing with previous measurements with equal names. Options are selfexplanatory and 'overwrite' is the default value.
            
        INPUT:
            - self: containing self.universe and self.measurements (if not loading a configuration file)

        OUPTUT:
            - self.results: dictionary containing arrays of the results keyed by the given name of the configuration.
        """

        if len(self.measurements) == 0:
            raise EmptyMeasurementsError
        else :
            for measurement in self.measurements:
                self.results[measurement["name"]] = []

        if last == -1: last = len(self.universe.trajectory) # so if no last frame given, last in trajectory is chosen.

        first = True
        for ts in tqdm(self.universe.trajectory[first-1:last:step], desc="Analysing", unit="frames"):
            print(first)
            for measurement in self.measurements:

                if measurement["type"] == "distance":

                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []

                    self.results[measurement["name"]].append(
                        distance(
                            measurement["sel"][0],
                            measurement["sel"][1],
                            measurement["options"]["type"]
                        )
                    )

                elif measurement["type"] == "dihedral":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []

                    self.results[measurement["name"]].append(
                        dihedral(
                            measurement["sel"][0],
                            measurement["sel"][1],
                            measurement["sel"][2],
                            measurement["sel"][3],
                            measurement["options"]["units"],
                            measurement["options"]["domain"]
                        )
                    )

                elif measurement["type"] == "angle":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        
                    self.results[measurement["name"]].append(
                        angle(
                            measurement["sel"][0],
                            measurement["sel"][1],
                            measurement["sel"][2],
                            measurement["options"]["units"],
                            measurement["options"]["domain"]
                        )
                    )

                elif measurement["type"] == "planar_angle":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        
                    self.results[measurement["name"]].append(
                        planar_angle(
                            plane_A=measurement["sel"][0].positions,
                            plane_B=measurement["sel"][1].positions,
                            units=measurement["options"]["units"],
                            domain=measurement["options"]["domain"],
                        )
                    )

                elif measurement["type"] == "contacts":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        
                    if measurement["mode"] == "selection":
                        self.results[measurement["name"]].append(contacts_selection(
                            sel= measurement["sel"][0],
                            sel_env= measurement["sel"][1],
                            interactions= measurement["options"]["interactions"],
                            out_format= measurement["options"]["out_format"],
                            measure_distances= measurement["options"]["measure_dists"],
                        ))

                    if measurement["mode"] == "protein":
                        self.results[measurement["name"]].append(contacts_protein(
                            sel= measurement["sel"][0],
                            sel_env= measurement["sel"][1],
                            interactions= measurement["options"]["interactions"],
                            out_format= measurement["options"]["out_format"],
                            measure_distances= measurement["options"]["measure_dists"],
                        ))

                    else :
                        pass

                elif measurement["type"] == "distWATbridge":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        
                    if (
                        "WAT" not in measurement["sel"][2].resnames
                        or "WAT" not in measurement["sel"][3].resnames
                    ) or (
                        len(
                            set(measurement["sel"][2].resids)
                            & set(measurement["sel"][3].resids)
                        )
                        == 0
                    ):  # No WAT in one or both selections' environment
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

                            if dist1_ <= measurement["sel"][4] and dist2_ <= measurement["sel"][5]:

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

                        selWAT = selection(
                                self.universe, int(closestWAT), sel_type="res_num"
                            )

                        self.results[measurement["name"]].append(
                            [closestWAT,
                            npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][0].positions),
                                    array(selWAT.positions),
                                    backend="OpenMP",
                                )
                            ),
                            npmin(
                                mdadist.distance_array(
                                    array(measurement["sel"][1].positions),
                                    array(selWAT.positions),
                                    backend="OpenMP",
                                )
                            )])

                        del selWAT, closestWAT, dist1_, dist1, dist2_, dist2

                elif measurement["type"] == "rmsd":
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        

                    if measurement["options"]["superposition"] == True:

                        self.results[measurement["name"]].append(
                            rms.rmsd(
                                measurement["sel"].positions
                                - measurement["sel"].center_of_mass(),
                                measurement["ref"],
                                center=True,
                                superposition=True,
                            )
                        )

                    elif measurement["options"]["superposition"] == False:

                        self.results[measurement["name"]].append(
                            rms.rmsd(
                                measurement["sel"].positions
                                - measurement["sel"].center_of_mass(),
                                measurement["ref"],
                            )
                        )

                elif measurement["type"] == 'pka':
                    if measurement["name"] in self.results and previous_measurements.lower() in ['ow', 'overwrite'] and first:
                        self.results[measurement["name"]] = []
                        
                    if first:
                        check_folder(measurement["options"]["pdb_folder"])

                    self.results[measurement["name"]].append(
                        pka(
                            sel_protein= measurement["sel"],
                            pka_ref= measurement["options"]["pka_ref"],
                            pdb_folder= measurement["options"]["pdb_folder"],
                            frame= str(ts.frame),
                            keep_pdb= measurement["options"]["keep_pdb"],
                            keep_pka= measurement["options"]["keep_pka"]
                        )
                    )

            if first:
                first = False

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
