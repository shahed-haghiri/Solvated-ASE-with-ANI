from ase.lattice.cubic import Diamond
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase import units
from ase import Atoms
import ase
import torchani
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.build import bulk, molecule
from ase.build.attach import attach_randomly, attach_randomly_and_broadcast
from ase.io.trajectory import Trajectory
from ase.io import Trajectory, write
import sys

# function that takes in a SMILES string and returns an ASE atoms object with coordinates from rdkit
# probably want to minimize it before using it
def smiles2atoms(smiles = 'C'):
    m = Chem.MolFromSmiles(smiles)
    m_Hs = Chem.AddHs(m)
    n_atoms = m_Hs.GetNumAtoms()
    AllChem.EmbedMolecule(m_Hs,randomSeed=0xf00d)
    coords3d = Chem.MolToMolBlock(m_Hs).split('\n')[4:4+n_atoms]
    coords = []
    atoms_string = []
    for i in coords3d:
        coords.append(tuple([float(x) for x in i.split()[:3]]))
        atoms_string.append(i.split()[3])
    atoms_string = ''.join(atoms_string)
    box = Atoms(atoms_string, positions=coords)
    return(box)

# runs minimization on given ase atoms object
def runminim(atoms_obj):
    opt = BFGS(atoms_obj)
    opt.run(fmax=0.001)
    print()

# adds chosen solvent randomly around the given atoms object n times 
# then returns the obj
def addsolvent(atoms_obj,solvent_obj,num_mol):
    for i in range(num_mol):
        atoms_obj = attach_randomly_and_broadcast(atoms_obj,solvent_obj,2)
    return(atoms_obj)

# function taken from torchANI ase tutorial
# https://aiqm.github.io/torchani/examples/ase_interface.html
# runs md for n steps and prints to log file and also makes a trajectory to be used by ASE
def runmd(atoms_obj, n_steps, temp, traj_name):
    def printenergy(a=atoms_obj):
        """Function to print the potential, kinetic and total energy."""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
            'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
        
    dyn = Langevin(atoms_obj, 1 * units.fs, temp * units.kB, 0.2,append_trajectory=True,trajectory=f'{traj_name}.traj',logfile=f'{traj_name}_log.txt')
    dyn.attach(printenergy, interval=50)

    print("Beginning dynamics...")
    printenergy()
    dyn.run(n_steps)


# function to convert the ASE trajectory into an xyz trajectory
def traj2xyz(traj_name):
    traj = Trajectory(f'{traj_name}.traj')
    write(f'{traj_name}_traj.xyz',traj)

def main():
    # input should be python3 solvated_ani_ase.py smiles solvent_smiles num_molecules temp steps
    smiles, solvent_smiles = sys.argv[1:3]
    num_molecules = int(sys.argv[3]) if len(sys.argv) >= 4 else 10
    temp = float(sys.argv[4]) if len(sys.argv) >= 5 else 298
    steps = int(sys.argv[5]) if len(sys.argv) >= 6 else 500 
    traj_name = f'{smiles}_{num_molecules}{solvent_smiles}_{int(temp)}K_{steps}_steps'
    
    box = smiles2atoms(smiles = smiles)
    solvent = smiles2atoms(smiles = solvent_smiles)
    
    # minimize both the solute and solvent before putting them together
    calculator = torchani.models.ANI1ccx().ase()
    box.set_calculator(calculator)
    runminim(box)
    solvent.set_calculator(calculator)
    runminim(solvent)

    # equilibrate the box
    box = addsolvent(box,solvent,num_molecules)
    box.set_calculator(calculator)
    runminim(box)
    
    # run md and convert to xyz trajectory
    runmd(box, n_steps = steps, temp = temp, traj_name = traj_name)
    traj2xyz(traj_name)

if __name__ == '__main__':
    main()
