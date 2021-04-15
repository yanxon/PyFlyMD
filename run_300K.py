import os
from ase import units
from ase.io import Trajectory
from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT

from pyxtal_ff import PyXtal_FF
from pyxtal_ff.calculator import PyXtalFFCalculator

from pyflymd import *

ase_db, dumpfolder = "Si.db", "vasp_Si"
if not os.path.exists(dumpfolder):
    os.makedirs(dumpfolder)
if os.path.isfile(ase_db):
    os.remove(ase_db)
kspacing = 0.2
steps = 10

struc = bulk("Si", 'diamond', a=5.468728, cubic=True)
mliap = "20-20-checkpoint.pth"
ff = PyXtal_FF(model={'system': ["Si"]}, logo=False)
ff.run(mode='predict', mliap=mliap)
calc = PyXtalFFCalculator(ff=ff)
struc.set_calculator(calc)

#struc.set_calculator(set_vasp('single', kspacing))
#energy, forces = dft_run(struc, path=dumpfolder, ncpu=16, max_time=10) # total energy

# MD
traj = Trajectory("Si.traj", 'w', struc)
MaxwellBoltzmannDistribution(struc, temperature_K=300)
dyn = NPT(struc, 0.1 * units.fs, externalstress=0., ttime=10*units.fs, pfactor=10*units.fs, 
          temperature_K=300)
dyn.attach(traj.write, interval=steps)

for i in range(100):
    dyn.run(steps=steps)
    PrintEnergy(struc)
    save_to_ASE_db(struc, path=ase_db)
traj.close()
