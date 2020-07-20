from dataclasses import dataclass

class QMInput:

    @dataclass
    class InputFiles:
        coordinates   : str
        topology      : str
        pdb           : str
        coords_chemsh : str
        active_atoms  : str


    @dataclass
    class CalculationFeatures:
        type           : str  = 'opt'
        optimiser      : str  = 'lbfgs'
        coordinates    : str  = 'hdlc'
        qmmm_embedding : str  = 'electrostatic'
        microiterative : bool = False
        inner_residues : list = None

    @dataclass
    class QM:
        qm_region                : str  = None
        software                 : str  = 'turbomole'
        hamiltonian              : str  = 'b3lyp'
        electronic_configuration : str  = 'closed'
        basis                    : str  = '6-31G(d,p)'
        metal                    : bool = False
        metal_basis              : str  = None
        charge                   : int  = 0
        multiplicity             : int  = 1
        extra                    : str  = None


    def create(self):
        """
            Function to create the input chemsh file from the class attributes
        """
        pass

    def print_input_files(self):
        print(self.InputFilescoordinates, self.InputFiles.topology, self.InputFiles.pdb, self.InputFiles.coords_chemsh, self.InputFiles.active_atoms)