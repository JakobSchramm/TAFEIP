# TAFEIP (Topological Aromaticity For Estimating Interface Properties)

Python Scripts for estimating the interface properties on a Cu(111) surface solely from the molecular structure using topological aromaticity (conjugated-circuit resonance energy & topological resonance energy) as a descriptor.
It is the result of two publications. If you use any of the scripts, please cite our work.

## Dependencies

- argparse
- ase
- ast
- itertools
- numpy
- pandas
- rdkit
- sympy

## Individual Scripts

- Estimate_Interface_Properties.py: Main script for estimating interface properties from molecular structure. Can be used directly from command line. Uses Generate_Connectivity.py, CCRE.py, and TRE.py.

- Generate_Connectivity.py: Script for generating connectivity (topological) matrix from .xyz-file, SMILES string, or matrix.
- Find_Rings.py: Script for finding all possible rings (circuits) in a structure (graph).
- Kekule.py: Script for generating all possible Kekule structures (chemical structure with explicit single and double bonds) for given connectivity (topological) matrix.
- CCRE.py: Script for calculating the conjugated-circuit resonance energy (CCRE). Uses Kekule.py and Find_Rings.py.
- TRE.py: Script for calculating the topological resonance energy (TRE). Uses Find_Rings.py.

## References

- J. Schramm, R. Tonner-Zech, Unveiling Correlations in Metal-Organic Interface Properties: A Computational Exploration of Alternant and Non-Alternant π-Electron Systems, ChemPlusChem 2025, 90 (6), e202400771, DOI: 10.1002/cplu.202400771.
- J. Schramm, R. Tonner-Zech, Molecular Aromaticity Determines Metal-Organic Interface Properties: Alternant and Non-Alternant π-Electron Systems on Cu(111), J. Phys. Chem. C 2025, submitted.
