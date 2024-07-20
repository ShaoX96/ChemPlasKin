# Introduction
paserPlasKin.py is part of ChemPlasKin project.

A python parser translating plasma kinetics input in ZDPlasKin format to YAML format.
Utilized to build a unified chemistry-plasma kinetic mechanism input.

Inline definition of variables will also be parsed, substituted and simplified.

Written by Xiao Shao @KAUST, 2023

# How to use
1. Open terminal
2. Ensure you have Python installed.
3. **Usage:**
```sh
python parsePlasKin.py --input <inputfile> [--output <outputfile>] [--markIndex <True/False>] [--quiet <True/False>]
```
4. **Example:**
```sh
python parsePlasKin.py --input "plasmaH2O2.inp" --output "parsedPlasKin.yaml" --markIndex True --quiet False
```
In a simple way:
```sh
python parsePlasKin.py --input "plasmaH2O2.inp"
```
6. If any of the following packages is missing in your python environment, install it:
- re
- yaml
- Counter
- sympy
- argparse
- logging

## What next:
Manually incorporate the parsed YAML content into the combustion mechanism files supported by Cantera.
Ensure the consistency of species names and existence of thermodynamic properties of plasma species.

# Note
Feel free to make modifications for your customized needs. If sticky errors are happening, try manually modify the output YAML file.