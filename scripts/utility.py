import paths, hashlib
import scm.plams as plams
import numpy as np


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


def get_sorted_dict_values(d):
    k, v = list(d.keys()), list(d.values())
    idx = np.argsort(list(k))
    return [v[i] for i in idx]



def hartree2kcalmol(v):
    return plams.Units.convert(v, 'Hartree', 'kcal/mol')