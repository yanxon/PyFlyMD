from ase.io import Trajectory
from ase.db import connect
from ase import units

def PrintEnergy(atoms):
    """Function to print the potential, kinetic and total energy"""
    epot = atoms.get_potential_energy() / len(atoms)
    ekin = atoms.get_kinetic_energy() / len(atoms)
    try:
        aleatoric = atoms.calc.results['aleatoric']/len(atoms)
        epistemic = atoms.calc.results['epistemic']/len(atoms)
    except:
        aleatoric = -1.
        epistemic = -1.
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV Alea = %.3feV Epis = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin, aleatoric, epistemic))

def save_to_ASE_db(atoms, path='ASE.db'):
    with connect(path) as db:
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()
        aleatoric = atoms.calc.results['aleatoric']
        epistemic = atoms.calc.results['epistemic']
        data = {'potential_energy': epot,
                'kinetic_energy': ekin,
                'total_energy': epot+ekin,
                'temperature': ekin / (1.5 * units.kB),
                'aleatoric': aleatoric,
                'epistemic': epistemic
                }
        db.write(atoms, data=data)

def print_trajectory(trajectory):
    data = Trajectory(trajectory)
    for atoms in data:
        PrintEnergy(atoms)
