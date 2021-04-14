import os
from ase.calculators.vasp import Vasp


def dft_run(struc, path, ncpu=1, max_time=3, clean=True):
    """
    perform dft calculation and get energy and forces
    """
    cmd = "mpirun -np " + str(ncpu) + " vasp_std"
    os.environ["VASP_COMMAND"] = "timeout " + str(max_time) + "m " + cmd
    cwd = os.getcwd()
    os.chdir(path)
    try:
        eng = struc.get_potential_energy()
        forces = struc.get_forces()
    except:
        eng = None
        forces = None
    if clean:
        os.system("rm POSCAR POTCAR INCAR OUTCAR")
    os.chdir(cwd)
    return eng, forces


def set_vasp(level='opt', kspacing=0.5):
    """
    INCAR parameters for VASP calculation
    """
    para0 = {'xc': 'pbe',
            'npar': 8,
            'kgamma': True,
            'lcharg': False,
            'lwave': False,
            'ibrion': 2,
            }
    if level == 'single':
        para = {'prec': 'accurate',
                'encut': 500,
                'ediff': 1e-4,
                'nsw': 0,
                }
    else:
        para = {'prec': 'accurate',
                'encut': 400,
                'isif': 3,
                'ediff': 1e-4,
                'nsw': 20,
                }
    dict_vasp = dict(para0, **para)
    return Vasp(kspacing=kspacing, **dict_vasp)
