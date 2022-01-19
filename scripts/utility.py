import paths
import scm.plams


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
	with open(path, 'r') as file:
		lines = file.readlines()
		flags = lines[1].split(',')
		mol.flags = [f.strip() for f in flags]

	return mol
