# CHANGELOG

## 0.3.0
- New description of the package both in setup.py and README.md
- Added the option to specify the first and last frames, and step to analyse in m.run_measure().
- New EmptyMeasurementsError, which raises when no measurement is given to run_measure().
- Defined errors/exceptions have been modified so Exception.__init__(self, 'error message') is used instead of print('error message') for better printing.
- Added the remove_measurement() function to Measurements() class to remove a measurement.
- New pka calculator for predicting the pKa for each tritable residue along the trajectory. The measurement can be added using the add_pKa function of the Measurements class.
- 