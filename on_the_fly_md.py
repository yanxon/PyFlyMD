import os
import shutil
from ase import units
from ase.io import Trajectory
from ase.build import bulk
from ase.db import connect
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase.md.langevin import Langevin

from pyxtal_ff import PyXtal_FF
from pyxtal_ff.calculator import PyXtalFFCalculator

from pyflymd import *

# Initial PyXtal_FF Training
TrainData = "Si_ase.db"
descriptor = {'type': 'SNAP',
              'Rc': 4.9,
              'weights': {'Si': 1.0},
              'parameters': {'lmax': 3},
              'ncpu': 1,
             }
model = {'system' : ['Si'],
         'hiddenlayers': [20, 20],
         'random_seed': 12345,
         'activation': ['Tanh', 'Tanh', 'Linear'],
         'optimizer': {'method': 'lbfgs'},
         'force_coefficient': 1.,
         'stress_coefficient': 1.,
         'alpha': 1e-6,
         'epoch': 500,
         }
ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode='train', TrainData=TrainData, TestData=None)
shutil.copyfile('Si-SNAP/20-20-checkpoint.pth', "potentials/Si-20-20-0K.pth")
shutil.rmtree('Si-SNAP')

Ts = [500, 750, 1000, 1250, 1500]
for i, T in enumerate(Ts):
    print("Temperature: ", T)
    ase_db = f"Si{T}.db"
    struc = bulk("Si", 'diamond', a=5.468728, cubic=True) * (2,2,2)
    struc.info['md_step'] = 0
    mliap = f'potentials/Si-20-20-0K.pth' if i == 0 else f'potentials/Si-20-20-{Ts[i-1]}K.pth'
    ff = PyXtal_FF(model={'system': ["Si"]}, logo=False)
    ff.run(mode='predict', mliap=mliap)
    calc = PyXtalFFCalculator(ff=ff)
    struc.set_calculator(calc)

    traj = Trajectory("Si.traj", "w", struc)
    MaxwellBoltzmannDistribution(struc, temperature_K=T)
    dyn = Langevin(struc, timestep=1 * units.fs, temperature_K=T, friction=0.08)
    dyn.attach(traj.write, interval=10)

    for _ in range(200):
        dyn.run(steps=10)
        PrintEnergy(struc)
        save_to_ASE_db(struc, path=ase_db)
    traj.close()
    
    with connect("Si_ase.db") as db:
        dft_structures = get_structures_with_high_aleatoric(f"Si{T}.db", size=1)
        for struc in dft_structures:
            struc.set_calculator(set_vasp('single', 0.2))
            eng, forces, stress = dft_run(struc, path = 'vasp_Si', ncpu=32, max_time=60)
            #print(eng)
            #print(type(eng))
            #print(forces)
            #print(stress)
            struc.set_calculator(None)
            stress = stress[[0,1,2,5,4,3]]
            data = {'energy': eng, 'force': forces, 'stress': stress, 'group': f'md{T}'}
            db.write(struc, data=data)

    model['restart'] = f'potentials/Si-20-20-0K.pth' if i == 0 else f'potentials/Si-20-20-{T}K.pth'
    ff = PyXtal_FF(descriptors=descriptor, model=model)
    ff.run(mode='train', TrainData=TrainData, TestData=None)
    shutil.copyfile('Si-SNAP/20-20-checkpoint.pth', f"potentials/Si-20-20-{T}K.pth")
    shutil.rmtree('Si-SNAP')
