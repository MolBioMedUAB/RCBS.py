from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()


setup(
    name="RCBS.py",
    version="0.3.0a",
    description="Python package for analysing MDs trajectories with MDAnalysis with complex analysing functions.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Miquel Canyelles Niño",
    author_mail="mcanyellesnino@gmail.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "pyyaml", "MDAnalysis", "tqdm"],
    keywords="biochemistry, simulations, MDAnalysis, molecular dynamics",
    url="https://github.com/dynamicsUAB/RCBS.py",
    download_url="https://github.com/dynamicsUAB/RCBS.py/archive/refs/tags/0.2.0.tar.gz",
)
