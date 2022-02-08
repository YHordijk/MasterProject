import os, sys
import periodictable as pt


def csv_to_dict(file):
	d = {}
	with open(file, 'r') as f:
		for row in f.readlines():
			row = row.strip('\n')
			row = row.split(', ')
			d[int(row[0])] = row[1:]

	return d

this_file_path = os.path.dirname(os.path.realpath(__file__))
RESOURCES_DIR = this_file_path + r'\data\resources'
if  not os.path.exists(RESOURCES_DIR):
	RESOURCES_DIR = this_file_path + r'\data\resources'

BASIS_SET_DIR = rf'{RESOURCES_DIR}\basissets'
XYZ_PATH = rf'{RESOURCES_DIR}\xyz'

#element parameters
MAX_PRINCIPAL_QUANTUM_NUMBER = csv_to_dict(rf'{RESOURCES_DIR}\elements\max_quantum_number.csv')
MAX_VALENCE = csv_to_dict(rf'{RESOURCES_DIR}\elements\max_valence.csv')
ATOM_COLOURS = csv_to_dict(rf'{RESOURCES_DIR}\elements\colours.csv')
ATOM_COLOURS = {i: [int(c[0]), int(c[1]), int(c[2])] for i, c in ATOM_COLOURS.items()}
IONISATION_ENERGIES = csv_to_dict(rf'{RESOURCES_DIR}\elements\ionisation_energies.csv')
COVALENT_RADII = {i: pt.elements[i].covalent_radius for i in range(96)}
ELECTRONEGATIVITIES = csv_to_dict(rf'{RESOURCES_DIR}\elements\electronegativities.csv')



#forcefield parameters
def load_uff_parameters(path):
	keys = ['Atom', 'r1', 'theta0', 'x1', 'D1', 'zeta', 'Z1', 'Vi', 'Uj', 'Xi', 'Hard', 'Radius']
	params = {k:[] for k in keys}
	with open(path, 'r') as prm:
		for l in prm.readlines():
			if l.startswith('param'):
				s = l.strip().split()[1:]
				[params[k].append(s[i]) for i, k in enumerate(keys)]

	return params


FF_UFF_PARAMETERS = load_uff_parameters(rf'{RESOURCES_DIR}\forcefields\UFF\UFF.prm')


def get_basis_file(name, extension='json'):
	f = rf'{BASIS_SET_DIR}\{name}.{extension}'
	assert(os.path.exists(f))
	return f

# print(get_basis_file('STO-4G'))