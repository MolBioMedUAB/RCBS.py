

class QMInput:
    def __init__(self):
        self.input_files          = self.InputFiles()
        self.calculation_features = self.CalculationFeatures()


    class InputFiles:
        def __init__(self):
            self.coordinates   = None
            self.topology      = None
            self.pdb           = None
            self.coords_chemsh = None
            self.active_atoms  = None


    class CalculationFeatures:
        def __init__(self):
            self.type           = 'opt'
            self.optimiser      = 'lbfgs'
            self.coordinates    = 'hdlc'
            self.qmmm_embedding = 'electrostatic'
            self.microiterative = False
            self.inner_residues = None

            self.qm = self.QM()

        class QM:
            def __init__(self):
                self.software                 = 'turbomole'
                self.qm_region                = None
                self.hamiltonian              = 'b3lyp'
                self.electronic_configuration = 'closed'
                self.basis                    = '6-31G(d,p)'
                self.metal                    = False
                self.metal_basis              = None
                self.charge                   = 0
                self.multiplicity             = 1
                self.extra                    = None


    def create():
        """
            Function to create the input chemsh file from the class attributes
        """
        pass