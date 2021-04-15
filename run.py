import os
import numpy as np
from ase import Atoms
from ase import units
from ase.io import Trajectory
from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT

from pyxtal_ff import PyXtal_FF
from pyxtal_ff.calculator import PyXtalFFCalculator

from pyflymd import *

dumpfolder = "vasp_Si"
if not os.path.exists(dumpfolder):
    os.makedirs(dumpfolder)
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
dyn = NPT(struc, 0.1 * units.fs, externalstress=0., ttime=5000*units.fs, pfactor=5000*units.fs, 
          temperature_K=300)
dyn.attach(traj.write, interval=steps)

for i in range(10):
    dyn.run(steps=steps)
    PrintEnergy(struc)
    save_to_ASE_db(struc, path='Si.db')

#traj.close()


#print(energy/len(struc))
#print(forces)
