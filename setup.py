from setuptools import setup, find_packages

with open("README.md", 'r') as f:
    long_description = f.read()


setup(
        name='RCBS.py',
        version='0.0',
        description='Python package useful for analysing MD trajectories and creating QM/MM models built on top of MDAnalysis. All MD formats compatible with MDAnalysis are compatible with RCBS, while QM/MM capacities are compatible with ChemShell. The aim of this package is to develop an easy-to-use way to analyse MD trajectories, while keeping a good performance.',
        long_description=long_description,
        author='Miquel Canyelles Niño',
        author_mail='mcanyellesnino@gmail.com',
        packages= find_packages(),
        include_package_data=True,
        install_requires=['numpy', 'pyyaml', 'MDAnalysis'])

