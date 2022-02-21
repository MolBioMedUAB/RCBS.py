from MDAnalysis import Universe
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.lib.distances as mdadist
from ..exceptions import NotEnoughAtomsSetectedError, NotSingleAtomSelectionError
from numpy import min as npmin
from numpy import array, matrix
from numpy.linalg import norm



class Measurements:
    def __init__(self, u):
        self.measurements = []
        #self.options      = []
        #self.type         = []
        self.universe     = u
        #self.trajectory   = u.trajectory
        self.results      = {}
        self.boolean      = {}


    def add_distance(self, name, sel1, sel2):
        """
        DESCRIPTION:
            This function outputs the minimum measured distance between the two input selections or coordinates or their combination.
        INPUT:
            - Two selections, which can contain more than one atom.
        OUTPUT:
            - Shorter distance between sel1 and sel2 (in ang).
        """

        self.measurements.append(
            {
                'name'    : name,
                'type'    : 'distance',
                'sel'     : [sel1, sel2],
                'options' : None,
 #               'measure' : None,
            }
            #npmin(mdadist.distance_array(array(sel1), array(sel2), backend='OpenMP'))
            )


    def add_dihedral(self, name, sel1, sel2, sel3, sel4, units='degree', domain=360):
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
            - Selection of four atoms in four different AtomGroups. They have to be input with the correct order.

        OUTPUT:
            - Dihedral angle between the input atoms.
        """

        for sel in (sel1, sel2, sel3, sel4):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError

        units = units.lower()
        domain = str(domain).lower()

        if units not in ('deg', 'degree', 'degrees', 'rad', 'radian', 'radians'):
            units = 'degree'

        if domain not in (180, '180', 360, '360', 'pi', '2pi'):
            domain = '360'

        self.measurements.append({
                'name'    : name,
                'type'    : 'dihedral',
                'sel'     : [sel1, sel2, sel3, sel4],
                'options' : {
                    'units'  : units,
                    'domain' : domain
                }
        })


    def add_angle(self, name, sel1, sel2, sel3, units='deg', domain=360):
        """
        DESCRIPTION:
            This functions measures the angle between 3 specified atoms and returns the value between 0 and 360 degrees.
            The input selections have to be single atoms.

        OPTIONS:
            - units: option for selecting the output units of the dihedral
                - degree
                - rad

            - domain: option for specifying the domain of the output measures
                - 180, pi: option for -180,180 domain
                - 360, 2pi: option for 0,360 domain. Default option

        INPUT:
            - Selection of four atoms in three different AtomGroups. They have to be input with the correct order.

        OUTPUT:
            - Angle between the input atoms.
        """


        for sel in (sel1, sel2, sel3):
            if len(sel) != 1:
                raise NotSingleAtomSelectionError

        units = units.lower()
        domain = str(domain).lower()

        if units not in ('deg', 'degree', 'degrees', 'rad', 'radian', 'radians'):
            units = 'degree'

        if domain not in (180, '180', 360, '360', 'pi', '2pi'):
            domain = '360'

        self.measurements.append({
                'name'    : name,
                'type'    : 'angle',
                'sel'     : [sel1, sel2, sel3],
                'options' : {
                    'units'  : units,
                    'domain' : domain
                }
        })



    def add_contacts(self, name, sel, sel_env=3, interactions='all', include_WAT=False):
        """
        DESCRIPTION:
            This function takes a Universe, a selection and a radius and returns the list of residues nearer than the specified radius.

        INPUT:
            - u            -> MDAnalysis Universe
            - sel          -> selection of central atoms. It has to be an AtomGroup
            - sel_env      -> radius (in ang)
            - interactions -> type of interactions to be considered (all, polar, nonpolar, donorHbond, none). Custom
                              interactions can be also analysed by passing a list of residues names.

        OUTPUT:
            - List of dictionaries containing the name and number of all interacting residues.
        """

        sel_env = self.universe.select_atoms('around %s group select' % sel_env, select=sel, updating=True)

        if isinstance(interactions, str):
            if interactions not in ('all', 'polar', 'nonpolar', 'donorHbond', 'none'):
                print("This type of interaction is not described. Available interactions are: 'all', 'polar', 'nonpolar', 'donorHbond' and 'none'.")
                return

            else :
                if interactions == 'all':
                    interactions = [
                        'ARG', 'HIS', 'HID', 'HIE', 'HIP', 'LYS',
                        'ASP', 'GLU',
                        'SER', 'THR', 'ASN', 'GLN',
                        'CYS', 'SEC', 'GLY', 'PRO',
                        'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL'
                    ]
                elif interactions == 'polar':
                    interactions = [
                        'ARG', 'HIS', 'HID', 'HIE', 'HIP', 'LYS',
                        'ASP', 'GLU',
                        'SER', 'THR', 'ASN', 'GLN',
                        'CYS', 'SEC',
                        'TYR'
                    ]

                elif interactions == 'nonpolar':
                    interactions = [
                        'CYS', 'SEC', 'GLY', 'PRO',
                        'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL'
                    ]

                elif interactions == 'donorHbond':
                    interactions = [
                        'ARG', 'HID', 'HIE', 'HIP', 'LYS',
                        'ASH', 'GLH',
                        'SER', 'THR', 'ASN', 'GLN',
                        'CYS', 'SEC', 'GLY', 'PRO',
                        'TYR'
                    ]

                elif interactions == 'none':
                    interactions = []

        elif isinstance(interactions, list):

            pass

        if include_WAT == True:
            interactions += ['WAT']


        self.measurements.append({
                'name'    : name,
                'type'    : 'contacts',
                'sel'     : [sel, sel_env],
                'options' : {
                    'interactions' : interactions,
                }
        })



    def config_saver(self, config_filename):
        """
        DESCRIPTION:
            Function for saving the configuration of a measurement run into a file (in json or yaml format).

        INPUT:
            - self (containing the measurements list of dictionaries) and the file name

        OUTPUT:
            - the dictionary saved in the file

        TODO:
            - [] Check what happens when selections are not AtomGroups but indices or names
        """

        if config_filename.split('.')[-1].lower() not in ('json', 'jsn', 'yaml', 'yml'):
            print("Output format file is not available. Use '.json' or '.yaml' extensions for saving in either format.")
            return

        #dict_ = []
        #for entry in range(len(self.measurements)):
        #    dict_.append(self.measurements[entry])

        config = []
        for l in self.measurements:
            config.append(l.copy())

        for measurement in config:

            sels = []
            for s in measurement['sel']:
                if isinstance(s, AtomGroup):
                    sels.append(s.indices)

                    for m in range(len(sels)):
                        sels[m] = list(sels[m])
                        for n in range(len(sels[m])):
                            sels[m][n] = int(sels[m][n] + 1)

                measurement['sel'] = sels

            measurement['sel_type'] = 'at_num'


        if config_filename.split('.')[-1].lower() in ('json', 'jsn'):
            from json import dump

            f = open(config_filename, 'w')
            dump(config, f)
            f.close()

            print('Configuration saved in', config_filename)

            del dump

        elif config_filename.split('.')[-1].lower() in ('yaml', 'yml'):
            from yaml import dump

            f = open(config_filename, 'w')
            dump(config, f)
            f.close()

            print('Configuration saved in', config_filename)

            del dump

        return config


    #def load_results():

    def run_measure(self, save_output=False, input_config=False):
        """

        """

        #results = {}
        for measurement in self.measurements:
            self.results[measurement['name']] = []

        for ts in self.universe.trajectory:

            for measurement in self.measurements:

                if measurement['type'] == 'distance':
                    self.results[measurement['name']].append(npmin(mdadist.distance_array(array(measurement['sel'][0].positions), array(measurement['sel'][1].positions), backend='OpenMP')))

                elif measurement['type'] == 'dihedral':

                    if measurement['options']['units'] in ('rad', 'radian', 'radians'):
                        if measurement['options']['domain'] in (180, '180', 'pi'):
                            self.results[measurement['name']].append(float(mdadist.calc_dihedrals(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, measurement['sel'][3].positions, backend='OpenMP')))
                        elif measurement['options']['domain'] in (360, '360', '2pi'):
                            from math import pi
                            self.results[measurement['name']].append((float(mdadist.calc_dihedrals(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, measurement['sel'][3].positions, backend='OpenMP')) + pi) % pi)

                    elif measurement['options']['units'] in ('deg', 'degree', 'degrees'):
                        from numpy import rad2deg
                        if measurement['options']['domain'] in (180, '180', 'pi'):
                            self.results[measurement['name']].append(float(rad2deg(mdadist.calc_dihedrals(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, measurement['sel'][3].positions, backend='OpenMP'))))
                        elif measurement['options']['domain'] in (360, '360', '2pi'):
                            from math import pi
                            self.results[measurement['name']].append((float(rad2deg(mdadist.calc_dihedrals(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, measurement['sel'][3].positions, backend='OpenMP'))) + 360) % 360)

                elif measurement['type'] == 'angle':

                    if measurement['options']['units'] in ('rad', 'radian', 'radians'):
                        if measurement['options']['domain'] in (180, '180', 'pi'):
                            self.results[measurement['name']].append(float(mdadist.calc_angles(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, backend='OpenMP')))
                        elif measurement['options']['domain'] in (360, '360', '2pi'):
                            from math import pi
                            self.results[measurement['name']].append((float(mdadist.calc_angles(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, backend='OpenMP')) + pi) % pi)

                    elif measurement['options']['units'] in ('deg', 'degree', 'degrees'):
                        from numpy import rad2deg

                        if measurement['options']['domain'] in (180, '180', 'pi'):
                            self.results[measurement['name']].append(float(rad2deg(mdadist.calc_angles(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, backend='OpenMP'))))
                        elif measurement['options']['domain'] in (360, '360', '2pi'):
                            from math import pi
                            self.results[measurement['name']].append((float(rad2deg(mdadist.calc_angles(measurement['sel'][0].positions, measurement['sel'][1].positions, measurement['sel'][2].positions, backend='OpenMP'))) + 360) % 360)



        if save_output != False:

            if save_output.split('.')[-1] in ('json', 'jsn'):
                from json import dump

                f = open(save_output, 'w')
                dump(self.results, f)
                f.close()

                print('Results saved in', save_output)

            elif save_output.split('.')[-1] in ('yaml', 'yml'):
                from yaml import dump

                f = open(save_output, 'w')
                dump(self.results, f)
                f.close()

                print('Results saved in', save_output)

        return self.results

    #def print_results_names(self):



    #    print()

    def run_boolean(self, *bool_configs, combine=True, save_output=False):
        """
        DESCRIPTION:
            Function for checking if

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
        """

        def bool_configs_checker(bool_configs):

            counter = 0
            for c in bool_configs:

                if 'measure_name' not in list(c.keys()):
                    c['measure_name'] = 'boolean_' + str(counter)

                if 'measure_type' not in list(c.keys()):
                    c['measure_type'] = 'dist'
                else:
                    if c['measure_type'].lower() not in ('dist', 'ang', 'dihe'):
                        c['measure_type'] = 'dist'

                if 'mode' not in list(c.keys()):
                    c['mode'] = 'lim'
                else :
                    if c['mode'].lower() not in ('lim', 'tol'):
                        c['mode'] = 'lim'

                if 'ref_val1' not in list(c.keys()):
                    while True:
                        try:
                            c['ref_val1'] = float(input('Input max (if lim mode) or central (if tol mode) value [float] for %s: ' % c['measure_name']))
                            break

                        except ValueError:
                            print('Input a number (int or float).')
                            continue

                if 'ref_val2' not in list(c.keys()):
                    while True:
                        try:
                            c['ref_val2'] = float(input('Input max (if lim mode) or central (if tol mode) value [float] for %s: ' % c['measure_name']))
                            break

                        except ValueError:
                            print('Input a number (int or float).')
                            continue

                counter += 1

            return bool_configs

        bool_configs = bool_configs_checker(bool_configs)


        for b in bool_configs:
            if b['measure_type'] in ('ang', 'dihe'):
                b['ref_val1'], b['ref_val2'] = list(map(lambda ang: ((ang + 360) % 360), [b['ref_val1'], b['ref_val2']]))

            if b['mode'] == 'lim':
                if b['ref_val1'] >= b['ref_val2']:
                    max_val = b['ref_val1']
                    min_val = b['ref_val2']
                elif b['ref_val1'] < b['ref_val2']:
                    max_val = b['ref_val2']
                    min_val = b['ref_val1']

            elif b['mode'] == 'tol':
                max_val = b['ref_val1'] + b['ref_val2']
                min_val = b['ref_val1'] - b['ref_val2']

            self.boolean[b['measure_name']] = []
            for r in self.results[b['measure_name']]:
                self.boolean[b['measure_name']].append(bool(r <= max_val and r > min_val))


        if combine == True:
            self.boolean['combination'] = []
            for f in range(len(self.boolean[list(self.boolean.keys())[0]])):

                boolean = True
                for key in list(self.boolean.keys())[:-1]:
                    boolean = bool(boolean and self.boolean[key][f])

                self.boolean['combination'].append(boolean)

        if save_output != False:

            if save_output.split('.')[-1] in ('json', 'jsn'):
                from json import dump

                f = open(save_output, 'w')
                dump(self.boolean, f)
                f.close()

                print('Boolean results saved in', save_output)

            elif save_output.split('.')[-1] in ('yaml', 'yml'):
                from yaml import dump

                f = open(save_output, 'w')
                dump(self.boolean, f)
                f.close()

                print('Boolean_results saved in', save_output)


        return self.boolean













