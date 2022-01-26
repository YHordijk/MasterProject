import paths, hashlib
import scm.plams as plams


def write_mol(mol, path, comment=None):
    with open(path, 'w+') as file:
        file.write(f'{len(mol.atoms)}\n')
        if comment is None:
        	comment = ', '.join(mol.flags)
        file.write(comment + '\n')
        for a in mol.atoms:
            file.write(f'{a.symbol}\t{a.coords[0]: .8f}\t{a.coords[1]: .8f}\t{a.coords[2]: .8f}\n')


def load_mol(path):
	mol = plams.Molecule(path)
	mol.path = path
	with open(path, 'r') as file:
		lines = file.readlines()
		flags = lines[1].split(',')
		mol.flags = [f.strip() for f in flags]

	return mol


def hash(reaction, stationary_point, flags):
	hashstr = reaction + stationary_point + ' '.join(list(sorted(flags)))
	return hashlib.sha256(hashstr.encode('utf-8')).hexdigest()


def hartree2kcalmol(v):
	return plams.Units.convert(v, 'Hartree', 'kcal/mol')