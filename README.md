# Solvated-ASE-with-ANI
Small repo for code that I made to run solvated systems from smiles string in ASE.
Made this small code to take in two smiles strings for a solute and solvent then add the solute around the solvent and run MD using ANI as a calculator.
Pretty much a quick and dirty way to sample solvated structures, also ASE did not have a nice way to solvate with any molecule so I used rdkit.

How to use:

run in command line with 5 arguments:
```
$python3 solvated_ase_ani.py [smiles] [solvent_smiles] [num_molecules] [temp] [steps]
```
- [smiles]: is the smiles string of the solute molecule you want
- [solvent_smiles]: is the smiles string of the solvent that you want 
## these first two are required or the script will not work ##
- [num_molecules]: integer for the number of solvent molecules you want to add, defaults to 10 if nothing is given
- [temp]: temperature that you want to run the MD at, in kelvin, defaults to 298 if nothing is given
- [steps]: number of steps that you want to run the MD for, defaults to 500 if nothing is given

### these next three are optional. but if you include one you have to include them all, all or none

so if you run something like:
```
$python3 solvated_ase_ani.py CC O 50 1000 100
```
it will generate a structure of one ethane in 50 waters, and run the MD at 1000K for 100 steps

Citations:
* TorchANI: https://aiqm.github.io/torchani/
* TorchANI ASE tutorial: https://aiqm.github.io/torchani/examples/ase_interface.html
* ASE: https://wiki.fysik.dtu.dk/ase/index.html
* rdkit: https://www.rdkit.org/docs/ 

need all of these installed before running the script
