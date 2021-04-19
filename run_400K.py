import os
from ase import units
from ase.io import Trajectory
from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase.md.langevin import Langevin

from pyxtal_ff import PyXtal_FF
from pyxtal_ff.calculator import PyXtalFFCalculator

from pyflymd import *

ase_db = "Si.db"
if os.path.isfile(ase_db):
    os.remove(ase_db)
steps = 1000
temp = 1000

struc = bulk("Si", 'diamond', a=5.468728, cubic=True) * (2,2,2)
struc.info['md_step'] = 0
mliap = "potentials/Si-20-20.pth"
ff = PyXtal_FF(model={'system': ["Si"]}, logo=False)
ff.run(mode='predict', mliap=mliap)
calc = PyXtalFFCalculator(ff=ff)
struc.set_calculator(calc)

def printenergy(atoms=struc):
    """Function to print the potential, kinetic and total energy"""
    step = atoms.info['md_step']
    epot = atoms.get_potential_energy() / len(atoms)
    ekin = atoms.get_kinetic_energy() / len(atoms)
    try:
        aleatoric = atoms.calc.results['aleatoric']/len(atoms)
        epistemic = atoms.calc.results['epistemic']/len(atoms)
    except:
        aleatoric = -1.
        epistemic = -1.
    print('%.4d Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV Alea = %.3feV Epis = %.3feV' % (step, epot, ekin, ekin / (1.5 * units.kB), epot + ekin, aleatoric, epistemic))
    atoms.info['md_step'] += 1

# MD
traj = Trajectory("Si.traj", 'w', struc)
MaxwellBoltzmannDistribution(struc, temperature_K=temp+50)
# NVT
dyn = Langevin(struc, timestep=1 * units.fs, temperature_K=temp, friction=0.08)
# NPT
#dyn = NPT(struc, 1. * units.fs, externalstress=0., ttime=1000*units.fs, pfactor=1000*units.fs, 
#          temperature_K=300)

dyn.attach(printenergy)
dyn.attach(traj.write, interval=10)
printenergy()
dyn.run(steps)

#for i in range(1000):
#    dyn.run(steps=steps)
#    PrintEnergy(struc, i)
#    save_to_ASE_db(struc, path=ase_db)
#traj.close()
