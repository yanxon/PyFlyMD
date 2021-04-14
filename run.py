import numpy as np
from ase import Atoms
from ase.build import bulk

from pyflymd import *

folder = "vasp_Si"

struc = bulk("Si", 'fcc', a=3.6, cubic=True)
struc.set_calculator(set_vasp('opt', 0.3))
eng, forces = dft_run(struc, path=folder, ncpu=4)
