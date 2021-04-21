from ase.db import connect

from pyflymd import *

path = 'Si.db'

Structures, Epot, Aleatoric, Epistemic, Temperatures = [], [], [], [], []
with connect(path) as db:
    for i, row, in enumerate(db.select()):
        struc = db.get_atoms(row.id)
        epot = row.data['potential_energy'] / len(struc)
        aleatoric = row.data["aleatoric"] / len(struc)
        epistemic = row.data["epistemic"] / len(struc)
        temperature = row.data["temperature"] / len(struc)

        Structures.append(struc)
        Epot.append(epot)
        Aleatoric.append(aleatoric)
        Epistemic.append(epistemic)
        Temperatures.append(temperature)

print(len(Structures))
import random

temp = 0.
while temp < 1000:
    num = 60 #random.randint(0, len(Structures))
    print("rand_num: ", num)
    struc = Structures[num]
    epot = Epot[num]
    temp = Temperatures[num]
    aleatoric = Aleatoric[num]
    epistemic = Epistemic[num]
print("Temperature: ", temp)

struc.set_calculator(set_vasp('single', 0.2))
energy, forces = dft_run(struc, path='vasp_Si', ncpu=32, max_time=60) # total energy

print("Temp:     ", temp)
print("Pred Eng: ", epot)
print("True Eng: ", energy/len(struc))
print("Diff Eng: ", abs(epot-energy/len(struc)))
print("Aleatoric: ", aleatoric)
print("Epistemic: ", epistemic)
