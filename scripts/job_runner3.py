import scm.plams as plams
import struct_generator2, paths, os, utility, time, sys, geometry
from utility import bohr2angstrom
b2a = bohr2angstrom(1)

join = os.path.join


if __name__ == '__main__':
    from job_runner import initialize, GO
    calc_dir = paths.calculations2

    all_calc_paths = []
    for R1 in ['Et', 'NMe2']:
        for R2 in ['Bz']:
            for cat in ['I2', 'ZnCl2', 'TiCl4', 'BF3', 'AlF3', 'SnCl4']:
                calc_path = initialize.initialize('achiral_radical_addition', {'R1':R1, 'R2':R2, 'Rcat':cat}, calc_dir=calc_dir)
                all_calc_paths.append(calc_path)


    for calc_path in all_calc_paths:
        for d in os.listdir(calc_path):
            if os.path.isdir(join(calc_path, d)):
                GO.GO(join(calc_path, d))
