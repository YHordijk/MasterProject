import paths, hashlib, os
import scm.plams as plams
import numpy as np


join = os.path.join

def write_mol(mol, path, comment=None):
    with open(path, 'w+') as file:
        file.write(f'{len(mol.atoms)}\n')
        if comment is None:
            comment = ', '.join(mol.flags)
        file.write(comment + '\n')
        for a in mol.atoms:
            file.write(f'{a.symbol}\t{a.coords[0]: .8f}\t{a.coords[1]: .8f}\t{a.coords[2]: .8f}\n')


def write_mol_list(l, path, comment=''):
    with open(path, 'w+') as file:
        file.write(f'{len(l)}\n')
        file.write(comment + '\n')
        for x in l:
            file.write(f'{x[0]:2}\t{x[1]:.8f}\t{x[2]:.8f}\t{x[3]:.8f}\n')


def load_mol(path):
    mol = plams.Molecule(path)
    mol.path = path
    with open(path, 'r') as file:
        lines = file.readlines()
        flags = lines[1].split(',')
        mol.flags = [f.strip() for f in flags]

    return mol


def hash(reaction, stationary_point, flags):
    # print(reaction, stationary_point, flags)
    hashstr = reaction + stationary_point + ' '.join(list(sorted(flags)))
    # print(hashstr)
    return hashlib.sha256(hashstr.encode('utf-8')).hexdigest()


def hash_from_info(info_path):
    with open(info_path, 'r') as info:
        lines = info.readlines()
        for line in lines:
            if line.startswith('unique='):
                hashstr = line.strip().split('=')[1]
                break
    return hashlib.sha256(hashstr.encode('utf-8')).hexdigest()


def hash2(hashstr):
    return hashlib.sha256(hashstr.encode('utf-8')).hexdigest()


def get_all_run_dirs(calc_dir):
    dirs = []
    for system in os.listdir(calc_dir):
        if not os.path.isdir(join(calc_dir, system)): continue

        for dir in os.listdir(join(calc_dir, system)):
            p = join(calc_dir,system,dir)
            dirs.append(p)
    return dirs


def hash_collision(h, calc_dir=paths.calculations):
    dirs = get_all_run_dirs(calc_dir)
    hashes = [hash_from_info(join(d, 'run.info')) for d in dirs]
    return h in hashes

def get_colliding_dirs(h, calc_dir=paths.calculations):
    dirs = get_all_run_dirs(calc_dir)
    dirs = [d for d in dirs if hash_from_info(join(d, 'run.info'))]
    return dirs



def get_sorted_dict_values(d):
    k, v = list(d.keys()), list(d.values())
    idx = np.argsort(list(k))
    return [v[i] for i in idx]



def hartree2kcalmol(v):
    return plams.Units.convert(v, 'Hartree', 'kcal/mol')

def bohr2angstrom(v):
    return plams.Units.convert(v, 'bohr', 'angstrom')


def print_table(header, data, tabs=0):
    #data is provided as columns
    assert len(header) == len(data)

    

