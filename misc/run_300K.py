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

ase_db, dumpVASP = "Si.db", "vasp_Si"
if os.path.isfile(ase_db):
    os.remove(ase_db)
#if not os.path.exists(dumpVASP):
#    os.makedirs(dumpVASP)
kspacing = 0.2
steps = 10
temp = 300

struc = bulk("Si", 'diamond', a=5.468728, cubic=True) * (3,3,3)
mliap = "potentials/Si-20-20.pth"
ff = PyXtal_FF(model={'system': ["Si"]}, logo=False)
ff.run(mode='predict', mliap=mliap)
calc = PyXtalFFCalculator(ff=ff)
struc.set_calculator(calc)

#struc.set_calculator(set_vasp('single', kspacing))
#energy, forces = dft_run(struc, path=dumpVASP, ncpu=16, max_time=10) # total energy

# MD
traj = Trajectory("Si.traj", 'w', struc)
MaxwellBoltzmannDistribution(struc, temperature_K=temp+50)
# NVT
dyn = Langevin(struc, timestep=0.15 * units.fs, temperature_K=temp, friction=0.004 * units.fs)
# NPT
#dyn = NPT(struc, 1. * units.fs, externalstress=0., ttime=1000*units.fs, pfactor=1000*units.fs, 
#          temperature_K=300)
dyn.attach(traj.write, interval=steps)

for i in range(1000):
    dyn.run(steps=steps)
    PrintEnergy(struc, i)
    save_to_ASE_db(struc, path=ase_db)
traj.close()
